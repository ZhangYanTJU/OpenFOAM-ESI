/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "constantFilmThermo.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(constantFilmThermo, 0);

addToRunTimeSelectionTable
(
    filmThermoModel,
    constantFilmThermo,
    dictionary
);


void constantFilmThermo::init(thermoData& td)
{
    if (coeffDict_.readIfPresent(td.name_, td.value_))
    {
        td.set_ = true;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

constantFilmThermo::constantFilmThermo
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    filmThermoModel(typeName, film, dict),
    name_(coeffDict_.lookup("specie")),
    rho0_("rho0"),
    mu0_("mu0"),
    sigma0_("sigma0"),
    Cp0_("Cp0"),
    kappa0_("kappa0"),
    D0_("D0"),
    hl0_("hl0"),
    pv0_("pv0"),
    W0_("W0"),
    Tb0_("Tb0")
{
    init(rho0_);
    init(mu0_);
    init(sigma0_);
    init(Cp0_);
    init(kappa0_);
    init(D0_);
    init(hl0_);
    init(pv0_);
    init(W0_);
    init(Tb0_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

constantFilmThermo::~constantFilmThermo()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const word& constantFilmThermo::name() const
{
    return name_;
}


scalar constantFilmThermo::rho
(
    const scalar p,
    const scalar T
) const
{
    if (!rho0_.set_)
    {
        coeffDict_.readEntry(rho0_.name_, rho0_.value_);
        rho0_.set_ = true;
    }

    return rho0_.value_;
}


scalar constantFilmThermo::mu
(
    const scalar p,
    const scalar T
) const
{
    if (!mu0_.set_)
    {
        coeffDict_.readEntry(mu0_.name_, mu0_.value_);
        mu0_.set_ = true;
    }

    return mu0_.value_;
}


scalar constantFilmThermo::sigma
(
    const scalar p,
    const scalar T
) const
{
    if (!sigma0_.set_)
    {
        coeffDict_.readEntry(sigma0_.name_, sigma0_.value_);
        sigma0_.set_ = true;
    }

    return sigma0_.value_;
}


scalar constantFilmThermo::Cp
(
    const scalar p,
    const scalar T
) const
{
    if (!Cp0_.set_)
    {
        coeffDict_.readEntry(Cp0_.name_, Cp0_.value_);
        Cp0_.set_ = true;
    }

    return Cp0_.value_;
}


scalar constantFilmThermo::kappa
(
    const scalar p,
    const scalar T
) const
{
    if (!kappa0_.set_)
    {
        coeffDict_.readEntry(kappa0_.name_, kappa0_.value_);
        kappa0_.set_ = true;
    }

    return kappa0_.value_;
}


scalar constantFilmThermo::D
(
    const scalar p,
    const scalar T
) const
{
    if (!D0_.set_)
    {
        coeffDict_.readEntry(D0_.name_, D0_.value_);
        D0_.set_ = true;
    }

    return D0_.value_;
}


scalar constantFilmThermo::hl
(
    const scalar p,
    const scalar T
) const
{
    if (!hl0_.set_)
    {
        coeffDict_.readEntry(hl0_.name_, hl0_.value_);
        hl0_.set_ = true;
    }

    return hl0_.value_;
}


scalar constantFilmThermo::pv
(
    const scalar p,
    const scalar T
) const
{
    if (!pv0_.set_)
    {
        coeffDict_.readEntry(pv0_.name_, pv0_.value_);
        pv0_.set_ = true;
    }

    return pv0_.value_;
}


scalar constantFilmThermo::W() const
{
    if (!W0_.set_)
    {
        coeffDict_.readEntry(W0_.name_, W0_.value_);
        W0_.set_ = true;
    }

    return W0_.value_;
}


scalar constantFilmThermo::Tb(const scalar p) const
{
    if (!Tb0_.set_)
    {
        coeffDict_.readEntry(Tb0_.name_, Tb0_.value_);
        Tb0_.set_ = true;
    }

    return Tb0_.value_;
}


tmp<volScalarField> constantFilmThermo::rho() const
{
    auto trho = volScalarField::New
    (
        IOobject::scopedName(type(), rho0_.name_),
        IOobject::NO_REGISTER,
        film().regionMesh(),
        dimensionedScalar(dimDensity, Zero),
        fvPatchFieldBase::extrapolatedCalculatedType()
    );

    trho.ref().primitiveFieldRef() = this->rho(0, 0);
    trho.ref().correctBoundaryConditions();

    return trho;
}


tmp<volScalarField> constantFilmThermo::mu() const
{
    auto tmu = volScalarField::New
    (
        IOobject::scopedName(type(), mu0_.name_),
        IOobject::NO_REGISTER,
        film().regionMesh(),
        dimensionedScalar(dimPressure*dimTime, Zero),
        extrapolatedCalculatedFvPatchScalarField::typeName
    );

    tmu.ref().primitiveFieldRef() = this->mu(0, 0);
    tmu.ref().correctBoundaryConditions();

    return tmu;
}


tmp<volScalarField> constantFilmThermo::sigma() const
{
    auto tsigma = volScalarField::New
    (
        IOobject::scopedName(type(), sigma0_.name_),
        IOobject::NO_REGISTER,
        film().regionMesh(),
        dimensionedScalar(dimMass/sqr(dimTime), Zero),
        extrapolatedCalculatedFvPatchScalarField::typeName
    );

    tsigma.ref().primitiveFieldRef() = this->sigma(0, 0);
    tsigma.ref().correctBoundaryConditions();

    return tsigma;
}


tmp<volScalarField> constantFilmThermo::Cp() const
{
    auto tCp = volScalarField::New
    (
        IOobject::scopedName(type(), Cp0_.name_),
        IOobject::NO_REGISTER,
        film().regionMesh(),
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, Zero),
        extrapolatedCalculatedFvPatchScalarField::typeName
    );

    tCp.ref().primitiveFieldRef() = this->Cp(0, 0);
    tCp.ref().correctBoundaryConditions();

    return tCp;
}


tmp<volScalarField> constantFilmThermo::kappa() const
{
    auto tkappa = volScalarField::New
    (
        IOobject::scopedName(type(), kappa0_.name_),
        IOobject::NO_REGISTER,
        film().regionMesh(),
        dimensionedScalar(dimPower/dimLength/dimTemperature, Zero),
        extrapolatedCalculatedFvPatchScalarField::typeName
    );

    tkappa.ref().primitiveFieldRef() = this->kappa(0, 0);
    tkappa.ref().correctBoundaryConditions();

    return tkappa;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
