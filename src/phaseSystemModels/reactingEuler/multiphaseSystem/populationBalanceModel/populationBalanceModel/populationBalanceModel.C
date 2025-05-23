/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2019 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "populationBalanceModel.H"
#include "coalescenceModel.H"
#include "breakupModel.H"
#include "binaryBreakupModel.H"
#include "driftModel.H"
#include "nucleationModel.H"
#include "phaseSystem.H"
#include "surfaceTensionModel.H"
#include "fvmDdt.H"
#include "fvcDdt.H"
#include "fvmSup.H"
#include "fvcSup.H"
#include "fvcDiv.H"
#include "phaseCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

void Foam::diameterModels::populationBalanceModel::registerVelocityGroups()
{
    forAll(fluid_.phases(), phasei)
    {
        if (isA<velocityGroup>(fluid_.phases()[phasei].dPtr()()))
        {
            const velocityGroup& velGroup =
                refCast<const velocityGroup>(fluid_.phases()[phasei].dPtr()());

            if (velGroup.popBalName() == this->name())
            {
                velocityGroups_.resize(velocityGroups_.size() + 1);

                velocityGroups_.set
                (
                    velocityGroups_.size() - 1,
                    &const_cast<velocityGroup&>(velGroup)
                );

                forAll(velGroup.sizeGroups(), i)
                {
                    this->registerSizeGroups
                    (
                        const_cast<sizeGroup&>(velGroup.sizeGroups()[i])
                    );
                }
            }
        }
    }
}


void Foam::diameterModels::populationBalanceModel::registerSizeGroups
(
    sizeGroup& group
)
{
    if
    (
        sizeGroups_.size() != 0
        &&
        group.x().value() <= sizeGroups_.last().x().value()
    )
    {
        FatalErrorInFunction
            << "Size groups must be entered according to their representative"
            << " size"
            << exit(FatalError);
    }

    sizeGroups_.resize(sizeGroups_.size() + 1);
    sizeGroups_.set(sizeGroups_.size() - 1, &group);

    // Grid generation over property space
    if (sizeGroups_.size() == 1)
    {
        // Set the first sizeGroup boundary to the representative value of
        // the first sizeGroup
        v_.append
        (
            new dimensionedScalar
            (
                "v",
                sizeGroups_.last().x()
            )
        );

        // Set the last sizeGroup boundary to the representative size of the
        // last sizeGroup
        v_.append
        (
            new dimensionedScalar
            (
                "v",
                sizeGroups_.last().x()
            )
        );
    }
    else
    {
        // Overwrite the next-to-last sizeGroup boundary to lie halfway between
        // the last two sizeGroups
        v_.last() =
            0.5
           *(
                sizeGroups_[sizeGroups_.size()-2].x()
              + sizeGroups_.last().x()
            );

        // Set the last sizeGroup boundary to the representative size of the
        // last sizeGroup
        v_.append
        (
            new dimensionedScalar
            (
                "v",
                sizeGroups_.last().x()
            )
        );
    }

    delta_.append(new PtrList<dimensionedScalar>());

    Su_.append
    (
        volScalarField::New
        (
            "Su",
            IOobject::NO_REGISTER,
            mesh_,
            dimensionedScalar(inv(dimTime), Zero)
        )
    );

    SuSp_.append
    (
        volScalarField::New
        (
            "SuSp",
            IOobject::NO_REGISTER,
            mesh_,
            dimensionedScalar(inv(dimTime), Zero)
        )
    );
}


void Foam::diameterModels::populationBalanceModel::createPhasePairs()
{
    forAll(velocityGroups_, i)
    {
        const phaseModel& phasei = velocityGroups_[i].phase();

        forAll(velocityGroups_, j)
        {
            const phaseModel& phasej = velocityGroups_[j].phase();

            if (&phasei != &phasej)
            {
                const phasePairKey key
                (
                    phasei.name(),
                    phasej.name(),
                    false
                );

                if (!phasePairs_.found(key))
                {
                    phasePairs_.insert
                    (
                        key,
                        autoPtr<phasePair>
                        (
                            new phasePair
                            (
                                phasei,
                                phasej
                            )
                        )
                    );
                }
            }
        }
    }
}


void Foam::diameterModels::populationBalanceModel::correct()
{
    calcDeltas();

    forAll(velocityGroups_, v)
    {
        velocityGroups_[v].preSolve();
    }

    forAll(coalescence_, model)
    {
        coalescence_[model].correct();
    }

    forAll(breakup_, model)
    {
        breakup_[model].correct();

        breakup_[model].dsdPtr()().correct();
    }

    forAll(binaryBreakup_, model)
    {
        binaryBreakup_[model].correct();
    }

    forAll(drift_, model)
    {
        drift_[model].correct();
    }

    forAll(nucleation_, model)
    {
        nucleation_[model].correct();
    }
}


void Foam::diameterModels::populationBalanceModel::
birthByCoalescence
(
    const label j,
    const label k
)
{
    const sizeGroup& fj = sizeGroups_[j];
    const sizeGroup& fk = sizeGroups_[k];

    dimensionedScalar Gamma;
    dimensionedScalar v = fj.x() + fk.x();

    for (label i = j; i < sizeGroups_.size(); i++)
    {
        // Calculate fraction for intra-interval events
        Gamma = gamma(i, v);

        if (Gamma.value() == 0) continue;

        const sizeGroup& fi = sizeGroups_[i];

        // Avoid double counting of events
        if (j == k)
        {
            Sui_ =
                0.5*fi.x()*coalescenceRate_()*Gamma
               *fj*fj.phase()/fj.x()
               *fk*fk.phase()/fk.x();
        }
        else
        {
            Sui_ =
                fi.x()*coalescenceRate_()*Gamma
               *fj*fj.phase()/fj.x()
               *fk*fk.phase()/fk.x();
        }

        Su_[i] += Sui_;

        const phasePairKey pairij
        (
            fi.phase().name(),
            fj.phase().name()
        );

        if (pDmdt_.found(pairij))
        {
            const scalar dmdtSign
            (
                Pair<word>::compare(pDmdt_.find(pairij).key(), pairij)
            );

            pDmdt_[pairij]->ref() += dmdtSign*fj.x()/v*Sui_*fi.phase().rho();
        }

        const phasePairKey pairik
        (
            fi.phase().name(),
            fk.phase().name()
        );

        if (pDmdt_.found(pairik))
        {
            const scalar dmdtSign
            (
                Pair<word>::compare(pDmdt_.find(pairik).key(), pairik)
            );

            pDmdt_[pairik]->ref() += dmdtSign*fk.x()/v*Sui_*fi.phase().rho();
        }
    }
}


void Foam::diameterModels::populationBalanceModel::
deathByCoalescence
(
    const label i,
    const label j
)
{
    const sizeGroup& fi = sizeGroups_[i];
    const sizeGroup& fj = sizeGroups_[j];

    SuSp_[i] += coalescenceRate_()*fi.phase()*fj*fj.phase()/fj.x();

    if (i != j)
    {
        SuSp_[j] += coalescenceRate_()*fj.phase()*fi*fi.phase()/fi.x();
    }
}


void Foam::diameterModels::populationBalanceModel::
birthByBreakup
(
    const label k,
    const label model
)
{
    const sizeGroup& fk = sizeGroups_[k];

    for (label i = 0; i <= k; i++)
    {
        const sizeGroup& fi = sizeGroups_[i];

        Sui_ =
            fi.x()*breakupRate_()*breakup_[model].dsdPtr()().nik(i, k)
           *fk*fk.phase()/fk.x();

        Su_[i] += Sui_;

        const phasePairKey pair
        (
            fi.phase().name(),
            fk.phase().name()
        );

        if (pDmdt_.found(pair))
        {
            const scalar dmdtSign
            (
                Pair<word>::compare(pDmdt_.find(pair).key(), pair)
            );

            pDmdt_[pair]->ref() += dmdtSign*Sui_*fi.phase().rho();
        }
    }
}


void Foam::diameterModels::populationBalanceModel::deathByBreakup(const label i)
{
    const sizeGroup& fi = sizeGroups_[i];

    SuSp_[i] += breakupRate_()*fi.phase();
}


void Foam::diameterModels::populationBalanceModel::calcDeltas()
{
    forAll(sizeGroups_, i)
    {
        if (delta_[i].empty())
        {
            for (label j = 0; j <= sizeGroups_.size() - 1; j++)
            {
                delta_[i].append
                (
                    new dimensionedScalar
                    (
                        "delta",
                        dimVolume,
                        v_[i+1].value() - v_[i].value()
                    )
                );

                if
                (
                    v_[i].value() < 0.5*sizeGroups_[j].x().value()
                 &&
                    0.5*sizeGroups_[j].x().value() < v_[i+1].value()
                )
                {
                    delta_[i][j] =  mag(0.5*sizeGroups_[j].x() - v_[i]);
                }
                else if (0.5*sizeGroups_[j].x().value() <= v_[i].value())
                {
                    delta_[i][j].value() = 0;
                }
            }
        }
    }
}


void Foam::diameterModels::populationBalanceModel::
birthByBinaryBreakup
(
    const label i,
    const label j
)
{
    const sizeGroup& fj = sizeGroups_[j];
    const sizeGroup& fi = sizeGroups_[i];

    Sui_ = fi.x()*binaryBreakupRate_()*delta_[i][j]*fj*fj.phase()/fj.x();

    Su_[i] += Sui_;

    const phasePairKey pairij
    (
        fi.phase().name(),
        fj.phase().name()
    );

    if (pDmdt_.found(pairij))
    {
        const scalar dmdtSign
        (
            Pair<word>::compare(pDmdt_.find(pairij).key(), pairij)
        );

        pDmdt_[pairij]->ref() += dmdtSign*Sui_*fi.phase().rho();
    }

    dimensionedScalar Gamma;
    dimensionedScalar v = fj.x() - fi.x();

    for (label k = 0; k <= j; k++)
    {
        // Calculate fraction for intra-interval events
        Gamma = gamma(k, v);

        if (Gamma.value() == 0) continue;

        const sizeGroup& fk = sizeGroups_[k];

        volScalarField& Suk = Sui_;

        Suk =
            sizeGroups_[k].x()*binaryBreakupRate_()*delta_[i][j]*Gamma
           *fj*fj.phase()/fj.x();

        Su_[k] += Suk;

        const phasePairKey pairkj
        (
            fk.phase().name(),
            fj.phase().name()
        );

        if (pDmdt_.found(pairkj))
        {
            const scalar dmdtSign
            (
                Pair<word>::compare
                (
                    pDmdt_.find(pairkj).key(),
                    pairkj
                )
            );

            pDmdt_[pairkj]->ref() += dmdtSign*Suk*fi.phase().rho();
        }
    }
}


void Foam::diameterModels::populationBalanceModel::
deathByBinaryBreakup
(
    const label j,
    const label i
)
{
    const volScalarField& alphai = sizeGroups_[i].phase();

    SuSp_[i] += alphai*binaryBreakupRate_()*delta_[j][i];
}


void Foam::diameterModels::populationBalanceModel::drift(const label i)
{
    const sizeGroup& fp = sizeGroups_[i];

    if (i == 0)
    {
        rx_() = pos(driftRate_())*sizeGroups_[i+1].x()/sizeGroups_[i].x();
    }
    else if (i == sizeGroups_.size() - 1)
    {
        rx_() = neg(driftRate_())*sizeGroups_[i-1].x()/sizeGroups_[i].x();
    }
    else
    {
        rx_() =
            pos(driftRate_())*sizeGroups_[i+1].x()/sizeGroups_[i].x()
          + neg(driftRate_())*sizeGroups_[i-1].x()/sizeGroups_[i].x();
    }

    SuSp_[i] +=
        (neg(1 - rx_()) + neg(rx_() - rx_()/(1 - rx_())))*driftRate_()
       *fp.phase()/((rx_() - 1)*fp.x());

    rx_() *= 0.0;
    rdx_() *= 0.0;

    if (i == sizeGroups_.size() - 2)
    {
        rx_() = pos(driftRate_())*sizeGroups_[i+1].x()/sizeGroups_[i].x();

        rdx_() =
            pos(driftRate_())
           *(sizeGroups_[i+1].x() - sizeGroups_[i].x())
           /(sizeGroups_[i].x() - sizeGroups_[i-1].x());
    }
    else if (i < sizeGroups_.size() - 2)
    {
        rx_() = pos(driftRate_())*sizeGroups_[i+2].x()/sizeGroups_[i+1].x();

        rdx_() =
            pos(driftRate_())
           *(sizeGroups_[i+2].x() - sizeGroups_[i+1].x())
           /(sizeGroups_[i+1].x() - sizeGroups_[i].x());
    }

    if (i == 1)
    {
        rx_() += neg(driftRate_())*sizeGroups_[i-1].x()/sizeGroups_[i].x();

        rdx_() +=
            neg(driftRate_())
           *(sizeGroups_[i].x() - sizeGroups_[i-1].x())
           /(sizeGroups_[i+1].x() - sizeGroups_[i].x());
    }
    else if (i > 1)
    {
        rx_() += neg(driftRate_())*sizeGroups_[i-2].x()/sizeGroups_[i-1].x();

        rdx_() +=
            neg(driftRate_())
           *(sizeGroups_[i-1].x() - sizeGroups_[i-2].x())
           /(sizeGroups_[i].x() - sizeGroups_[i-1].x());
    }

    if (i != sizeGroups_.size() - 1)
    {
        const sizeGroup& fe = sizeGroups_[i+1];
        volScalarField& Sue = Sui_;

        Sue =
            pos(driftRate_())*driftRate_()*rdx_()
           *fp*fp.phase()/fp.x()
           /(rx_() - 1);

        Su_[i+1] += Sue;

        const phasePairKey pairij
        (
            fp.phase().name(),
            fe.phase().name()
        );

        if (pDmdt_.found(pairij))
        {
            const scalar dmdtSign
            (
                Pair<word>::compare(pDmdt_.find(pairij).key(), pairij)
            );

            pDmdt_[pairij]->ref() -= dmdtSign*Sue*fp.phase().rho();
        }
    }

    if (i != 0)
    {
        const sizeGroup& fw = sizeGroups_[i-1];
        volScalarField& Suw = Sui_;

        Suw =
            neg(driftRate_())*driftRate_()*rdx_()
           *fp*fp.phase()/fp.x()
           /(rx_() - 1);

        Su_[i-1] += Suw;

        const phasePairKey pairih
        (
            fp.phase().name(),
            fw.phase().name()
        );

        if (pDmdt_.found(pairih))
        {
            const scalar dmdtSign
            (
                Pair<word>::compare(pDmdt_.find(pairih).key(), pairih)
            );

            pDmdt_[pairih]->ref() -= dmdtSign*Suw*fp.phase().rho();
        }
    }
}


void Foam::diameterModels::populationBalanceModel::nucleation(const label i)
{
    Su_[i] += sizeGroups_[i].x()*nucleationRate_();
}


void Foam::diameterModels::populationBalanceModel::sources()
{
    forAll(sizeGroups_, i)
    {
        Su_[i] *= 0.0;
        SuSp_[i] *= 0.0;
    }

    forAllConstIter
    (
        phasePairTable,
        phasePairs(),
        phasePairIter
    )
    {
        pDmdt_(phasePairIter())->ref() *= 0.0;
    }

    // Since the calculation of the rates is computationally expensive,
    // they are calculated once for each sizeGroup pair and inserted into source
    // terms as required
    forAll(sizeGroups_, i)
    {
        const sizeGroup& fi = sizeGroups_[i];

        if (coalescence_.size() != 0)
        {
            for (label j = 0; j <= i; j++)
            {
                const sizeGroup& fj = sizeGroups_[j];

                if (fi.x() + fj.x() > sizeGroups_.last().x()) break;

                coalescenceRate_() *= 0.0;

                forAll(coalescence_, model)
                {
                    coalescence_[model].addToCoalescenceRate
                    (
                        coalescenceRate_(),
                        i,
                        j
                    );
                }

                birthByCoalescence(i, j);

                deathByCoalescence(i, j);
            }
        }

        if (breakup_.size() != 0)
        {
            forAll(breakup_, model)
            {
                breakup_[model].setBreakupRate(breakupRate_(), i);

                birthByBreakup(i, model);

                deathByBreakup(i);
            }
        }

        if (binaryBreakup_.size() != 0)
        {
            label j = 0;

            while (delta_[j][i].value() != 0)
            {
                binaryBreakupRate_() *= 0.0;

                forAll(binaryBreakup_, model)
                {
                    binaryBreakup_[model].addToBinaryBreakupRate
                    (
                        binaryBreakupRate_(),
                        j,
                        i
                    );
                }

                birthByBinaryBreakup(j, i);

                deathByBinaryBreakup(j, i);

                j++;
            }
        }

        if (drift_.size() != 0)
        {
            driftRate_() *= 0.0;

            forAll(drift_, model)
            {
                drift_[model].addToDriftRate(driftRate_(), i);
            }

            drift(i);
        }

        if (nucleation_.size() != 0)
        {
            nucleationRate_() *= 0.0;

            forAll(nucleation_, model)
            {
                nucleation_[model].addToNucleationRate(nucleationRate_(), i);
            }

            nucleation(i);
        }
    }
}


void Foam::diameterModels::populationBalanceModel::dmdt()
{
    forAll(velocityGroups_, v)
    {
        velocityGroup& velGroup = velocityGroups_[v];

        velGroup.dmdtRef() *= 0.0;

        forAll(sizeGroups_, i)
        {
            if (&sizeGroups_[i].phase() == &velGroup.phase())
            {
                sizeGroup& fi = sizeGroups_[i];

                velGroup.dmdtRef() += fi.phase().rho()*(Su_[i] - SuSp_[i]*fi);
            }
        }
    }
}


void Foam::diameterModels::populationBalanceModel::calcAlphas()
{
    alphas_() *= 0.0;

    forAll(velocityGroups_, v)
    {
        const phaseModel& phase = velocityGroups_[v].phase();

        alphas_() += max(phase, phase.residualAlpha());
    }
}


Foam::tmp<Foam::volScalarField>
Foam::diameterModels::populationBalanceModel::calcDsm()
{
    tmp<volScalarField> tInvDsm
    (
        volScalarField::New
        (
            "invDsm",
            mesh_,
            dimensionedScalar("zero", inv(dimLength), Zero)
        )
    );

    volScalarField& invDsm = tInvDsm.ref();

    forAll(velocityGroups_, v)
    {
        const phaseModel& phase = velocityGroups_[v].phase();

        invDsm += max(phase, phase.residualAlpha())/(phase.d()*alphas_());
    }

    return 1.0/tInvDsm;
}


void Foam::diameterModels::populationBalanceModel::calcVelocity()
{
    U_() *= 0.0;

    forAll(velocityGroups_, v)
    {
        const phaseModel& phase = velocityGroups_[v].phase();

        U_() += max(phase, phase.residualAlpha())*phase.U()/alphas_();
    }
}

bool Foam::diameterModels::populationBalanceModel::updateSources()
{
    const bool result = sourceUpdateCounter_ % sourceUpdateInterval() == 0;

    ++ sourceUpdateCounter_;

    return result;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::populationBalanceModel::populationBalanceModel
(
    const phaseSystem& fluid,
    const word& name,
    HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash>& pDmdt
)
:
    regIOobject
    (
        IOobject
        (
            name,
            fluid.time().constant(),
            fluid.mesh()
        )
    ),
    fluid_(fluid),
    pDmdt_(pDmdt),
    mesh_(fluid_.mesh()),
    name_(name),
    dict_
    (
        fluid_.subDict("populationBalanceCoeffs").subDict(name_)
    ),
    pimple_(mesh_.lookupObject<pimpleControl>("solutionControl")),
    continuousPhase_
    (
        mesh_.lookupObject<phaseModel>
        (
            IOobject::groupName("alpha", dict_.get<word>("continuousPhase"))
        )
    ),
    velocityGroups_(),
    sizeGroups_(),
    v_(),
    delta_(),
    Su_(),
    SuSp_(),
    Sui_
    (
        IOobject
        (
            "Sui",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("zero", inv(dimTime), Zero)
    ),
    coalescence_
    (
        dict_.lookup("coalescenceModels"),
        coalescenceModel::iNew(*this)
    ),
    coalescenceRate_(),
    breakup_
    (
        dict_.lookup("breakupModels"),
        breakupModel::iNew(*this)
    ),
    breakupRate_(),
    binaryBreakup_
    (
        dict_.lookup("binaryBreakupModels"),
        binaryBreakupModel::iNew(*this)
    ),
    binaryBreakupRate_(),
    drift_
    (
        dict_.lookup("driftModels"),
        driftModel::iNew(*this)
    ),
    driftRate_(),
    rx_(),
    rdx_(),
    nucleation_
    (
        dict_.lookup("nucleationModels"),
        nucleationModel::iNew(*this)
    ),
    nucleationRate_(),
    alphas_(),
    dsm_(),
    U_(),
    sourceUpdateCounter_
    (
        (mesh_.time().timeIndex()*nCorr())%sourceUpdateInterval()
    )
{
    this->registerVelocityGroups();

    this->createPhasePairs();

    if (sizeGroups_.size() < 3)
    {
        FatalErrorInFunction
            << "The populationBalance " << name_
            << " requires a minimum number of three sizeGroups to be specified."
            << exit(FatalError);
    }


    if (coalescence_.size() != 0)
    {
        coalescenceRate_.reset
        (
            new volScalarField
            (
                mesh_.newIOobject("coalescenceRate"),
                mesh_,
                dimensionedScalar(dimVolume/dimTime, Zero)
            )
        );
    }

    if (breakup_.size() != 0)
    {
        breakupRate_.reset
        (
            new volScalarField
            (
                mesh_.newIOobject("breakupRate"),
                mesh_,
                dimensionedScalar(inv(dimTime), Zero)
            )
        );
    }

    if (binaryBreakup_.size() != 0)
    {
        binaryBreakupRate_.reset
        (
            new volScalarField
            (
                mesh_.newIOobject("binaryBreakupRate"),
                mesh_,
                dimensionedScalar(inv(dimVolume*dimTime), Zero)
            )
        );
    }

    if (drift_.size() != 0)
    {
        driftRate_.reset
        (
            new volScalarField
            (
                mesh_.newIOobject("driftRate"),
                mesh_,
                dimensionedScalar(dimVolume/dimTime, Zero)
            )
        );

        rx_.reset
        (
            new volScalarField
            (
                mesh_.newIOobject("rx"),
                mesh_,
                dimensionedScalar(dimless, Zero)
            )
        );

        rdx_.reset
        (
            new volScalarField
            (
                mesh_.newIOobject("rdx"),
                mesh_,
                dimensionedScalar(dimless, Zero)
            )
        );
    }

    if (nucleation_.size() != 0)
    {
        nucleationRate_.reset
        (
            new volScalarField
            (
                mesh_.newIOobject("nucleationRate"),
                mesh_,
                dimensionedScalar(inv(dimTime*dimVolume), Zero)
            )
        );
    }

    if (velocityGroups_.size() > 1)
    {
        alphas_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("alpha", this->name()),
                    fluid_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    IOobject::REGISTER
                ),
                mesh_,
                dimensionedScalar(dimless, Zero)
            )
        );

        dsm_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("d", this->name()),
                    fluid_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    IOobject::REGISTER
                ),
                mesh_,
                dimensionedScalar(dimLength, Zero)
            )
        );

        U_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    IOobject::groupName("U", this->name()),
                    fluid_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    IOobject::REGISTER
                ),
                mesh_,
                dimensionedVector(dimVelocity, Zero)
            )
        );
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::populationBalanceModel::~populationBalanceModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::diameterModels::populationBalanceModel>
Foam::diameterModels::populationBalanceModel::clone() const
{
    NotImplemented;
    return nullptr;
}


bool Foam::diameterModels::populationBalanceModel::writeData(Ostream& os) const
{
    return os.good();
}


const Foam::dimensionedScalar
Foam::diameterModels::populationBalanceModel::gamma
(
    const label i,
    const dimensionedScalar& v
) const
{
    dimensionedScalar lowerBoundary(v);
    dimensionedScalar upperBoundary(v);
    const dimensionedScalar& xi = sizeGroups_[i].x();

    if (i == 0)
    {
       lowerBoundary = xi;
    }
    else
    {
       lowerBoundary = sizeGroups_[i-1].x();
    }

    if (i == sizeGroups_.size() - 1)
    {
        upperBoundary = xi;
    }
    else
    {
        upperBoundary = sizeGroups_[i+1].x();
    }

    if (v < lowerBoundary || v > upperBoundary)
    {
        return 0;
    }
    else if (v.value() <= xi.value())
    {
        return (v - lowerBoundary)/(xi - lowerBoundary);
    }
    else
    {
        return (upperBoundary - v)/(upperBoundary - xi);
    }
}


const Foam::tmp<Foam::volScalarField>
Foam::diameterModels::populationBalanceModel::sigmaWithContinuousPhase
(
    const phaseModel& dispersedPhase
) const
{
    return
        fluid_.lookupSubModel<reactingMultiphaseEuler::surfaceTensionModel>
        (
            phasePair(dispersedPhase, continuousPhase_)
        ).sigma();
}


const Foam::phaseCompressibleTurbulenceModel&
Foam::diameterModels::populationBalanceModel::continuousTurbulence() const
{
    return
        mesh_.lookupObject<phaseCompressibleTurbulenceModel>
        (
            IOobject::groupName
            (
                turbulenceModel::propertiesName,
                continuousPhase_.name()
            )
        );
}


void Foam::diameterModels::populationBalanceModel::solve()
{
    const dictionary& solutionControls = mesh_.solverDict(name_);
    const bool solveOnFinalIterOnly =
        solutionControls.getOrDefault("solveOnFinalIterOnly", false);

    if (!solveOnFinalIterOnly || pimple_.finalIter())
    {
        const label nCorr = this->nCorr();
        const scalar tolerance = solutionControls.get<scalar>("tolerance");

        if (nCorr > 0)
        {
            correct();
        }

        int iCorr = 0;
        scalar maxInitialResidual = 1;

        while (++iCorr <= nCorr && maxInitialResidual > tolerance)
        {
            Info<< "populationBalance "
                << this->name()
                << ": Iteration "
                << iCorr
                << endl;

            if (updateSources())
            {
                sources();
            }

            dmdt();

            maxInitialResidual = 0;

            forAll(sizeGroups_, i)
            {
                sizeGroup& fi = sizeGroups_[i];
                const phaseModel& phase = fi.phase();
                const volScalarField& alpha = phase;
                const dimensionedScalar& residualAlpha = phase.residualAlpha();
                const tmp<volScalarField> trho(phase.thermo().rho());
                const volScalarField& rho = trho();

                fvScalarMatrix sizeGroupEqn
                (
                    fvm::ddt(alpha, rho, fi)
                  + fi.VelocityGroup().mvConvection()->fvmDiv
                    (
                        phase.alphaRhoPhi(),
                        fi
                    )
                  - fvm::Sp
                    (
                        fvc::ddt(alpha, rho) + fvc::div(phase.alphaRhoPhi())
                      - fi.VelocityGroup().dmdt(),
                        fi
                    )
                  ==
                    fvc::Su(Su_[i]*rho, fi)
                  - fvm::SuSp(SuSp_[i]*rho, fi)
                  + fvc::ddt(residualAlpha*rho, fi)
                  - fvm::ddt(residualAlpha*rho, fi)
                );

                sizeGroupEqn.relax();

                maxInitialResidual = max
                (
                    sizeGroupEqn.solve().initialResidual(),
                    maxInitialResidual
                );
            }
        }

        if (nCorr > 0)
        {
            forAll(velocityGroups_, i)
            {
                velocityGroups_[i].postSolve();
            }
        }

        if (velocityGroups_.size() > 1)
        {
            calcAlphas();
            dsm_() = calcDsm();
            calcVelocity();
        }

        volScalarField fAlpha0
        (
            sizeGroups_.first()*sizeGroups_.first().phase()
        );

        volScalarField fAlphaN
        (
            sizeGroups_.last()*sizeGroups_.last().phase()
        );

        Info<< this->name() << " sizeGroup phase fraction first, last = "
            << fAlpha0.weightedAverage(this->mesh().V()).value()
            << ' ' << fAlphaN.weightedAverage(this->mesh().V()).value()
            << endl;
    }
}

// ************************************************************************* //
