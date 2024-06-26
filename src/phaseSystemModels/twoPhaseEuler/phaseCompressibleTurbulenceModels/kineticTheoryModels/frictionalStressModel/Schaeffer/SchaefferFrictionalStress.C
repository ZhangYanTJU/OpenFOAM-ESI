/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "SchaefferFrictionalStress.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace frictionalStressModels
{
    defineTypeNameAndDebug(Schaeffer, 0);

    addToRunTimeSelectionTable
    (
        frictionalStressModel,
        Schaeffer,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::Schaeffer::Schaeffer
(
    const dictionary& dict
)
:
    frictionalStressModel(dict),
    coeffDict_(dict.optionalSubDict(typeName + "Coeffs")),
    phi_("phi", dimless, coeffDict_)
{
    phi_ *= degToRad();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::Schaeffer::~Schaeffer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::Schaeffer::
frictionalPressure
(
    const phaseModel& phase,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax
) const
{
    const volScalarField& alpha = phase;

    return
        dimensionedScalar("1e24", dimensionSet(1, -1, -2, 0, 0), 1e24)
       *pow(Foam::max(alpha - alphaMinFriction, scalar(0)), 10.0);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::Schaeffer::
frictionalPressurePrime
(
    const phaseModel& phase,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax
) const
{
    const volScalarField& alpha = phase;

    return
        dimensionedScalar("1e25", dimensionSet(1, -1, -2, 0, 0), 1e25)
       *pow(Foam::max(alpha - alphaMinFriction, scalar(0)), 9.0);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::Schaeffer::nu
(
    const phaseModel& phase,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax,
    const volScalarField& pf,
    const volSymmTensorField& D
) const
{
    const volScalarField& alpha = phase;

    auto tnu = volScalarField::New
    (
        IOobject::scopedName("Schaeffer", "nu"),
        IOobject::NO_REGISTER,
        phase.mesh(),
        dimensionedScalar(dimensionSet(0, 2, -1, 0, 0), Zero)
    );
    auto& nuf = tnu.ref();

    forAll(D, celli)
    {
        if (alpha[celli] > alphaMinFriction.value())
        {
            nuf[celli] =
                0.5*pf[celli]*sin(phi_.value())
               /(
                    sqrt((1.0/3.0)*sqr(tr(D[celli])) - invariantII(D[celli]))
                  + SMALL
                );
        }
    }

    const fvPatchList& patches = phase.mesh().boundary();
    const tmp<volVectorField> tU(phase.U());
    const volVectorField& U = tU();

    volScalarField::Boundary& nufBf = nuf.boundaryFieldRef();

    forAll(patches, patchi)
    {
        if (!patches[patchi].coupled())
        {
            nufBf[patchi] =
                (
                    pf.boundaryField()[patchi]*sin(phi_.value())
                   /(
                        mag(U.boundaryField()[patchi].snGrad())
                      + SMALL
                    )
                );
        }
    }

    // Correct coupled BCs
    nuf.correctBoundaryConditions();

    return tnu;
}


bool Foam::kineticTheoryModels::frictionalStressModels::Schaeffer::read()
{
    coeffDict_ <<= dict_.optionalSubDict(typeName + "Coeffs");

    phi_.read(coeffDict_);
    phi_ *= degToRad();

    return true;
}


// ************************************************************************* //
