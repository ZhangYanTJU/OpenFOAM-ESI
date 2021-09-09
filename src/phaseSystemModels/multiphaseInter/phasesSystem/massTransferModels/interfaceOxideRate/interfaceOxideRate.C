/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "interfaceOxideRate.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::meltingEvaporationModels::interfaceOxideRate<Thermo, OtherThermo>
::interfaceOxideRate
(
    const dictionary& dict,
    const phasePair& pair
)
:
    InterfaceCompositionModel<Thermo, OtherThermo>(dict, pair),
    C_
    (
        dimensionedScalar
        (
            dimDensity/dimTime,
            dict.getCheck<scalar>("C", scalarMinMax::ge(0))
        )
    ),
    Tliquidus_
    (
        dimensionedScalar
        (
            dimTemperature,
            dict.getCheck<scalar>("Tliquidus", scalarMinMax::ge(0))
        )
    ),
    Tsolidus_
    (
        dimensionedScalar
        (
            dimTemperature,
            dict.getCheck<scalar>("Tsolidus", scalarMinMax::ge(0))
        )
    ),
    oxideCrit_
    (
        dimensionedScalar
        (
            dimDensity,
            dict.getCheck<scalar>("oxideCrit", scalarMinMax::ge(0))
        )
    ),
    mDotOxide_
    (
        IOobject
        (
            "mDotOxide",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimDensity/dimTime, Zero)
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::interfaceOxideRate<Thermo, OtherThermo>::Kexp
(
    const volScalarField& T
)
{
    const volScalarField& from = this->pair().from();
    const volScalarField& to = this->pair().to();

    // (CSC:Eq. 2)
    tmp<volScalarField> Salpha = from*(1 - from);

    // (CSC:Eq. 5)
    tmp<volScalarField> Soxide =
        max((oxideCrit_.value() - to)/oxideCrit_.value(), scalar(0));

    // (CSC:Eq. 4)
    tmp<volScalarField> ST =
        Foam::exp(1 - 1/max((T - Tsolidus_)/(Tliquidus_ - Tsolidus_), 1e-6));

    // (CSC:Eq. 6)
    mDotOxide_ = C_*Salpha*Soxide*ST;

    return tmp<volScalarField>::New(mDotOxide_);
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::interfaceOxideRate<Thermo, OtherThermo>::KSp
(
    label variable,
    const volScalarField& refValue
)
{
    return nullptr;
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::interfaceOxideRate<Thermo, OtherThermo>::KSu
(
    label variable,
    const volScalarField& refValue
)
{
    return nullptr;
}


// ************************************************************************* //
