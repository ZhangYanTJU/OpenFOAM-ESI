/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "solidChemistryModel.H"
#include "reactingMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class SolidThermo>
Foam::solidChemistryModel<CompType, SolidThermo>::solidChemistryModel
(
    typename CompType::reactionThermo& thermo
)
:
    CompType(thermo),
    ODESystem(),
    Ys_(this->solidThermo().composition().Y()),
    reactions_
    (
        dynamic_cast<const reactingMixture<SolidThermo>&>
        (
            this->solidThermo()
        )
    ),
    solidThermo_
    (
        dynamic_cast<const reactingMixture<SolidThermo>&>
        (
            this->solidThermo()
        ).speciesData()
    ),
    nSolids_(Ys_.size()),
    nReaction_(reactions_.size()),
    RRs_(nSolids_),
    reactingCells_(this->mesh().nCells(), true)
{
    // create the fields for the chemistry sources
    forAll(RRs_, fieldi)
    {
        RRs_.set
        (
            fieldi,
            new volScalarField::Internal
            (
                IOobject
                (
                    "RRs." + Ys_[fieldi].name(),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimMass/dimVolume/dimTime, Zero)
            )
        );
   }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class SolidThermo>
Foam::solidChemistryModel<CompType, SolidThermo>::
~solidChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class SolidThermo>
Foam::scalar Foam::solidChemistryModel<CompType, SolidThermo>::solve
(
    const scalarField& deltaT
)
{
    NotImplemented;
    return 0;
}


template<class CompType, class SolidThermo>
Foam::tmp<Foam::volScalarField>
Foam::solidChemistryModel<CompType, SolidThermo>::tc() const
{
    NotImplemented;
    return volScalarField::null();
}


template<class CompType, class SolidThermo>
Foam::tmp<Foam::volScalarField>
Foam::solidChemistryModel<CompType, SolidThermo>::Qdot() const
{
    auto tQdot = volScalarField::New
    (
        "Qdot",
        IOobject::NO_REGISTER,
        this->mesh_,
        dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
    );
    scalarField& Qdot = tQdot.ref();

    if (this->chemistry_)
    {
        forAll(Ys_, i)
        {
            forAll(Qdot, celli)
            {
                scalar hf = solidThermo_[i].Hc();
                Qdot[celli] -= hf*RRs_[i][celli];
            }
        }
    }

    return tQdot;
}


template<class CompType, class SolidThermo>
void Foam::solidChemistryModel<CompType, SolidThermo>::setCellReacting
(
    const label celli,
    const bool active
)
{
    reactingCells_[celli] = active;
}

// ************************************************************************* //
