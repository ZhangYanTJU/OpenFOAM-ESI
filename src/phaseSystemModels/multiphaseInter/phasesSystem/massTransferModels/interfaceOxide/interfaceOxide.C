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

#include "interfaceOxide.H"
#include "constants.H"
#include "cutCellIso.H"
#include "volPointInterpolation.H"
#include "wallPolyPatch.H"
#include "fvcSmooth.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::meltingEvaporationModels::interfaceOxide<Thermo, OtherThermo>
::interfaceOxide
(
    const dictionary& dict,
    const phasePair& pair
)
:
    InterfaceCompositionModel<Thermo, OtherThermo>(dict, pair),
    C_("C", dimDensity/dimTime, dict),
    Tactivate_("Tactivate", dimTemperature, dict),
    alphaCrit_("alphaCrit", dimless, dict),

    mDotc_
    (
        IOobject
        (
            "mDotc",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimDensity/dimTime, Zero)
    )
{
    if (sign(C_.value()) < 0)
    {

    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::interfaceOxide<Thermo, OtherThermo>
::Kexp(const volScalarField& T)
{
    const volScalarField ST(Foam::exp(1 - 1/max(T/Tactivate_, 1e-6)));

    const volScalarField& from = this->pair().from();
    const volScalarField& to = this->pair().to();

    const volScalarField boundedAlphaSqr
    (
        from*(1-from)*max((alphaCrit_.value() - to)/alphaCrit_.value(), scalar(0))
    );

    mDotc_ = C_*boundedAlphaSqr*ST;

    return tmp<volScalarField>(new volScalarField(mDotc_));
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::interfaceOxide<Thermo, OtherThermo>::KSp
(
    label variable,
    const volScalarField& refValue
)
{
    return tmp<volScalarField> ();
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::interfaceOxide<Thermo, OtherThermo>::KSu
(
    label variable,
    const volScalarField& refValue
)
{
    return tmp<volScalarField> ();
}


template<class Thermo, class OtherThermo>
const Foam::dimensionedScalar&
Foam::meltingEvaporationModels::interfaceOxide<Thermo, OtherThermo>
::Tactivate() const
{
    return Tactivate_;
}


template<class Thermo, class OtherThermo>
bool
Foam::meltingEvaporationModels::
interfaceOxide<Thermo, OtherThermo>::includeDivU()
{
    return true;
}

// ************************************************************************* //
