/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "IsothermalPhaseModel.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::IsothermalPhaseModel<BasePhaseModel>::IsothermalPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, index)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::IsothermalPhaseModel<BasePhaseModel>::correctThermo()
{
    BasePhaseModel::correctThermo();

    // Correct the thermo, but make sure that the temperature remains the same
    auto TCopy = volScalarField::New
    (
        IOobject::scopedName(this->thermo().T().name(), "Copy"),
        IOobject::NO_REGISTER,
        this->thermo().T()
    );
    this->thermo_->he() = this->thermo().he(this->thermo().p(), TCopy());
    this->thermo_->correct();
    this->thermo_->T() = TCopy;
}


template<class BasePhaseModel>
bool Foam::IsothermalPhaseModel<BasePhaseModel>::isothermal() const
{
    return true;
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::IsothermalPhaseModel<BasePhaseModel>::heEqn()
{
    FatalErrorInFunction
        << "Cannot construct an energy equation for an isothermal phase"
        << exit(FatalError);

    return nullptr;
}


// ************************************************************************* //
