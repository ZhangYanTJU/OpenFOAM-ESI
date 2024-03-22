/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "grass.H"
#include "basicThermo.H"
#include "fluidThermo.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace evapotranspirationHeatTransferModels
{
    defineTypeNameAndDebug(grass, 0);
    addToRunTimeSelectionTable
    (
        evapotranspirationHeatTransferModel,
        grass,
        dictionary
    );
}
}


const Foam::Enum
<
    Foam::evapotranspirationHeatTransferModels::grass::soilHeatFluxType
>
Foam::evapotranspirationHeatTransferModels::grass::soilHeatFluxTypeNames
({
    {
        soilHeatFluxType::PROPORTIONAL_TO_SOLAR_RADIATION,
        "proportionalToSolarRadiation"
    },
    { soilHeatFluxType::BOUNDARY, "boundary" }
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::evapotranspirationHeatTransferModels::grass::E(const labelList& cells)
const
{
    const scalar q = this->q();
    const scalar G = this->G(cells);
    const scalar Delta = this->Delta();
    const scalar D = this->D();
    const scalar ra = this->ra();
    const scalar rs = this->rs();
    const scalar gamma = this->gamma();

    // (BSG:Eq. 11)
    const scalar E =
        (Delta*(q - G) + rho_*Cp_*D/ra)/(Delta + gamma*(1 + rs/ra));

    return tmp<scalarField>::New(cells.size(), E);
}


Foam::scalar Foam::evapotranspirationHeatTransferModels::grass::Delta() const
{
    // (BSG:Eq. 15), (APR:Eq. 13) - note 0.6108 -> 610.8 due to kPa -> Pa
    const scalar Ta = Tref_ + 237.3;

    return 4098*(610.8*exp(17.27*Tref_/Ta))/sqr(Ta);
}


Foam::scalar Foam::evapotranspirationHeatTransferModels::grass::D() const
{
    // (BSG:Eq. 16)
    const scalar RHref = RHptr_->value(mesh().time().timeOutputValue());

    return (1 - RHref)*pSat();
}


Foam::scalar Foam::evapotranspirationHeatTransferModels::grass::pSat() const
{
    // (BSG:Eq. 17; Arden Buck equation)
    const scalar p1 = (1.0007 + 3.46e-8*pAtm_)*611.21;

    return p1*exp((18.678 - Tref_/234.5)*Tref_/(Tref_ + 257.14));
}


Foam::scalar Foam::evapotranspirationHeatTransferModels::grass::gamma() const
{
    // (APR:Eq. 8)
    return Cp_*pAtm_/(epsilon_*lambda());
}


Foam::scalar Foam::evapotranspirationHeatTransferModels::grass::lambda() const
{
    // (BSG:Eq. 12)
    const scalar T1 = Tref_ + 273.15;
    const scalar T2 = Tref_ + 239.24;

    return 1.91846e6*sqr(T1/T2);
}


Foam::scalar Foam::evapotranspirationHeatTransferModels::grass::ra() const
{
    // (APR:Eq. 4), (BSG:Eq. 14)
    const scalar log1 = log((zRefU_ - d_*h_)/(zom_*h_));
    const scalar log2 = log((zRefH_ - d_*h_)/(zoh_*h_));

    return log1*log2/(sqr(kappa_)*uRef_);
}


Foam::scalar Foam::evapotranspirationHeatTransferModels::grass::rs() const
{
    // (BSG:Eq. 13), (APR:Eq. 5; reduced form in p. 4-5)
    return ri_/(scalar(12)*h_);
}


Foam::scalar Foam::evapotranspirationHeatTransferModels::grass::S
(
    const labelList& cells
) const
{
    // Mark fvOption cells within mesh
    bitSet isZoneCell(mesh().nCells());
    isZoneCell.set(cells);

    scalar S = 0;

    // Select cells next to specified patches
    // Sum patch area that is covered by fvOption cells
    for (const label patchi : patchSet_)
    {
        const scalarField& s = mesh().magSf().boundaryField()[patchi];

        const polyPatch& pp = mesh().boundaryMesh()[patchi];
        const labelList& faceCells = pp.faceCells();

        forAll(faceCells, i)
        {
            const bool isCovered = isZoneCell[faceCells[i]];

            if (isCovered)
            {
                S += s[i];
            }
        }
    }
    reduce(S, sumOp<scalar>());

    if (mag(S) < SMALL)
    {
        FatalErrorInFunction
            << "Area coverage of grass cannot be zero. "
            << "Check whether cellZone has any boundary faces."
            << exit(FatalError);
    }

    return S;
}


Foam::scalar Foam::evapotranspirationHeatTransferModels::grass::G
(
    const labelList& cells
) const
{
    switch (soilHeatFluxMethod_)
    {
        case soilHeatFluxType::PROPORTIONAL_TO_SOLAR_RADIATION:
        {
            return Csoil_*q();
        }

        case soilHeatFluxType::BOUNDARY:
        {
            // Retrieve heat flux through patches
            tmp<FieldField<Field, scalar>> tqBf = this->qBf();
            const FieldField<Field, scalar>& qBf = tqBf.cref();

            // Mark fvOption cells within mesh
            bitSet isZoneCell(mesh().nCells());
            isZoneCell.set(cells);

            scalar G = 0;

            for (const label patchi : patchSet_)
            {
                const scalarField& s = mesh().magSf().boundaryField()[patchi];

                const scalarField& qfld = qBf[patchi];

                const polyPatch& pp = mesh().boundaryMesh()[patchi];
                const labelList& faceCells = pp.faceCells();

                forAll(faceCells, i)
                {
                    const bool isCovered = isZoneCell[faceCells[i]];

                    if (isCovered)
                    {
                        G += qfld[i]*s[i];
                    }
                }
            }
            reduce(G, sumOp<scalar>());

            return G/area_;
        }
    }

    return -1;
}


Foam::tmp<Foam::FieldField<Foam::Field, Foam::scalar>>
Foam::evapotranspirationHeatTransferModels::grass::qBf() const
{
    const auto& T = mesh().lookupObject<volScalarField>(TName_);
    const volScalarField::Boundary& Tbf = T.boundaryField();

    auto tq = tmp<FieldField<Field, scalar>>::New(Tbf.size());
    auto& q = tq.ref();

    forAll(q, patchi)
    {
        q.set(patchi, new Field<scalar>(Tbf[patchi].size(), Zero));
    }

    typedef compressible::turbulenceModel cmpTurbModel;

    if (mesh().foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const auto& turb =
            mesh().lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

        // Note: calling he(p,T) instead of he()
        const volScalarField he(turb.transport().he(turb.transport().p(), T));
        const volScalarField::Boundary& hebf = he.boundaryField();

        const volScalarField alphaEff(turb.alphaEff());
        const volScalarField::Boundary& alphaEffbf = alphaEff.boundaryField();

        for (const label patchi : patchSet_)
        {
            q[patchi] = alphaEffbf[patchi]*hebf[patchi].snGrad();
        }
    }
    else if (mesh().foundObject<fluidThermo>(fluidThermo::dictName))
    {
        const auto& thermo =
            mesh().lookupObject<fluidThermo>(fluidThermo::dictName);

        // Note: calling he(p,T) instead of he()
        const volScalarField he(thermo.he(thermo.p(), T));
        const volScalarField::Boundary& hebf = he.boundaryField();

        const volScalarField& alpha(thermo.alpha());
        const volScalarField::Boundary& alphabf = alpha.boundaryField();

        for (const label patchi : patchSet_)
        {
            q[patchi] = alphabf[patchi]*hebf[patchi].snGrad();
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find a valid thermo model to evaluate q. " << nl
            << "Database contents are: " << mesh().objectRegistry::sortedToc()
            << exit(FatalError);
    }

    // No radiative heat flux contribution is present

    return tq;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::evapotranspirationHeatTransferModels::grass::grass
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    evapotranspirationHeatTransferModel(dict, mesh),
    soilHeatFluxMethod_
    (
        soilHeatFluxTypeNames.getOrDefault
        (
            "soilHeatFluxMethod",
            dict,
            soilHeatFluxType::PROPORTIONAL_TO_SOLAR_RADIATION
        )
    ),
    Tptr_(nullptr),
    RHptr_(nullptr),
    area_(-1),
    Csoil_(),
    rho_(),
    Cp_(),
    epsilon_(),
    h_(),
    kappa_(),
    uRef_(),
    zRefU_(),
    zRefH_(),
    zom_(),
    zoh_(),
    d_(),
    ri_(),
    pAtm_(),
    Tref_(0),
    TName_(),
    patchSet_()
{
    Info<< "    Activating evapotranspiration heat transfer model: "
        << typeName << endl;

    read(dict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::evapotranspirationHeatTransferModels::grass::Q
(
    const labelList& cells
) const
{
    if (area_ == -1)
    {
        area_ = S(cells);
    }

    Tref_ = Tptr_->value(mesh().time().timeOutputValue());

    return E(cells)*area_;
}


bool Foam::evapotranspirationHeatTransferModels::grass::read
(
    const dictionary& dict
)
{
    if (!evapotranspirationHeatTransferModel::read(dict))
    {
        return false;
    }

    Tptr_.reset(Function1<scalar>::New("Tref", dict, &mesh()));
    RHptr_.reset(Function1<scalar>::New("RHref", dict, &mesh()));

    area_ = -1;

    const auto range = scalarMinMax::ge(SMALL);

    Csoil_ = dict.getCheckOrDefault<scalar>("Csoil", 0.1, range);
    rho_ = dict.getCheckOrDefault<scalar>("rho", 1.225, range);
    Cp_ = dict.getCheckOrDefault<scalar>("Cp", 1013.0, range);
    epsilon_ = dict.getCheckOrDefault<scalar>("epsilon", 0.622, range);
    h_ = dict.getCheckOrDefault<scalar>("h", 0.1, range);
    kappa_ = dict.getCheckOrDefault<scalar>("kappa", 0.41, range);
    uRef_ = dict.getCheckOrDefault<scalar>("uRef", 2, range);
    zRefU_ = dict.getCheckOrDefault<scalar>("zRefU", 10, range);
    zRefH_ = dict.getCheckOrDefault<scalar>("zRefH", 10, range);
    zom_ = dict.getCheckOrDefault<scalar>("zom", 0.123, range);
    zoh_ = dict.getCheckOrDefault<scalar>("zoh", 0.0123, range);
    d_ = dict.getCheckOrDefault<scalar>("d", 2.0/3.0, range);
    ri_ = dict.getCheckOrDefault<scalar>("ri", 100, range);
    pAtm_ = dict.getCheckOrDefault<scalar>("pAtm", 101.325, range);
    TName_ = dict.getOrDefault<word>("T", "T");

    if (soilHeatFluxMethod_ == soilHeatFluxType::BOUNDARY)
    {
        patchSet_ =
            mesh().boundaryMesh().patchSet(dict.get<wordRes>("patches"));
    }

    return true;
}


// ************************************************************************* //
