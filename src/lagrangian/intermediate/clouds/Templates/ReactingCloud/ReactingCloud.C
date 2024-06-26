/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "ReactingCloud.H"
#include "CompositionModel.H"
#include "PhaseChangeModel.H"

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class CloudType>
void Foam::ReactingCloud<CloudType>::setModels()
{
    compositionModel_.reset
    (
        CompositionModel<ReactingCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    phaseChangeModel_.reset
    (
        PhaseChangeModel<ReactingCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
}


template<class CloudType>
void Foam::ReactingCloud<CloudType>::checkSuppliedComposition
(
    const scalarField& YSupplied,
    const scalarField& Y,
    const word& YName
)
{
    if (YSupplied.size() != Y.size())
    {
        FatalErrorInFunction
            << YName << " supplied, but size is not compatible with "
            << "parcel composition: " << nl << "    "
            << YName << "(" << YSupplied.size() << ") vs required composition "
            << YName << "(" << Y.size() << ")" << nl
            << abort(FatalError);
    }
}


template<class CloudType>
void Foam::ReactingCloud<CloudType>::cloudReset(ReactingCloud<CloudType>& c)
{
    CloudType::cloudReset(c);

    compositionModel_.reset(c.compositionModel_.ptr());
    phaseChangeModel_.reset(c.phaseChangeModel_.ptr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ReactingCloud<CloudType>::ReactingCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const SLGThermo& thermo,
    bool readFields
)
:
    CloudType(cloudName, rho, U, g, thermo, false),
    reactingCloud(),
    cloudCopyPtr_(nullptr),
    constProps_(this->particleProperties()),
    compositionModel_(nullptr),
    phaseChangeModel_(nullptr),
    rhoTrans_(thermo.carrier().species().size())
{
    if (this->solution().active())
    {
        setModels();

        if (readFields)
        {
            parcelType::readFields(*this, this->composition());
            this->deleteLostParticles();
        }
    }

    // Set storage for mass source fields and initialise to zero
    forAll(rhoTrans_, i)
    {
        const word fieldName
        (
            IOobject::scopedName
            (
                this->name(),
                "rhoTrans_" + thermo.carrier().species()[i]
            )
        );

        rhoTrans_.set
        (
            i,
            new volScalarField::Internal
            (
                IOobject
                (
                    fieldName,
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE,
                    IOobject::REGISTER
                ),
                this->mesh(),
                dimensionedScalar(dimMass, Zero)
            )
        );
    }

    if (this->solution().resetSourcesOnStartup())
    {
        resetSourceTerms();
    }
}


template<class CloudType>
Foam::ReactingCloud<CloudType>::ReactingCloud
(
    ReactingCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    reactingCloud(),
    cloudCopyPtr_(nullptr),
    constProps_(c.constProps_),
    compositionModel_(c.compositionModel_->clone()),
    phaseChangeModel_(c.phaseChangeModel_->clone()),
    rhoTrans_(c.rhoTrans_.size())
{
    forAll(c.rhoTrans_, i)
    {
        const word fieldName
        (
            IOobject::scopedName
            (
                this->name(),
                "rhoTrans_" + this->thermo().carrier().species()[i]
            )
        );

        rhoTrans_.set
        (
            i,
            new volScalarField::Internal
            (
                IOobject
                (
                    fieldName,
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                c.rhoTrans_[i]
            )
        );
    }
}


template<class CloudType>
Foam::ReactingCloud<CloudType>::ReactingCloud
(
    const fvMesh& mesh,
    const word& name,
    const ReactingCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    reactingCloud(),
    cloudCopyPtr_(nullptr),
    constProps_(),
    compositionModel_(c.compositionModel_->clone()),
    phaseChangeModel_(nullptr),
    rhoTrans_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ReactingCloud<CloudType>::setParcelThermoProperties
(
    parcelType& parcel,
    const scalar lagrangianDt
)
{
    CloudType::setParcelThermoProperties(parcel, lagrangianDt);

    parcel.Y() = composition().YMixture0();
}


template<class CloudType>
void Foam::ReactingCloud<CloudType>::checkParcelProperties
(
    parcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    CloudType::checkParcelProperties(parcel, lagrangianDt, fullyDescribed);

    if (fullyDescribed)
    {
        checkSuppliedComposition
        (
            parcel.Y(),
            composition().YMixture0(),
            "YMixture"
        );
    }

    // derived information - store initial mass
    parcel.mass0() = parcel.mass();
}


template<class CloudType>
void Foam::ReactingCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<ReactingCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::ReactingCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::ReactingCloud<CloudType>::resetSourceTerms()
{
    CloudType::resetSourceTerms();
    forAll(rhoTrans_, i)
    {
        rhoTrans_[i].field() = 0.0;
    }
}


template<class CloudType>
void Foam::ReactingCloud<CloudType>::relaxSources
(
    const ReactingCloud<CloudType>& cloudOldTime
)
{
    CloudType::relaxSources(cloudOldTime);

    typedef volScalarField::Internal dsfType;

    forAll(rhoTrans_, fieldi)
    {
        dsfType& rhoT = rhoTrans_[fieldi];
        const dsfType& rhoT0 = cloudOldTime.rhoTrans()[fieldi];
        this->relax(rhoT, rhoT0, "rho");
    }
}


template<class CloudType>
void Foam::ReactingCloud<CloudType>::scaleSources()
{
    CloudType::scaleSources();

    typedef volScalarField::Internal dsfType;

    forAll(rhoTrans_, fieldi)
    {
        dsfType& rhoT = rhoTrans_[fieldi];
        this->scale(rhoT, "rho");
    }
}


template<class CloudType>
void Foam::ReactingCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::trackingData td(*this);

        this->solve(*this, td);
    }
}


template<class CloudType>
void Foam::ReactingCloud<CloudType>::autoMap(const mapPolyMesh& mapper)
{
    Cloud<parcelType>::autoMap(mapper);

    this->updateMesh();
}


template<class CloudType>
void Foam::ReactingCloud<CloudType>::info()
{
    CloudType::info();

    this->phaseChange().info();
}


template<class CloudType>
void Foam::ReactingCloud<CloudType>::writeFields() const
{
    if (compositionModel_)
    {
        CloudType::particleType::writeFields(*this, this->composition());
    }
}


template<class CloudType>
void Foam::ReactingCloud<CloudType>::writeObjects(objectRegistry& obr) const
{
    CloudType::particleType::writeObjects(*this, this->composition(), obr);
}


// ************************************************************************* //
