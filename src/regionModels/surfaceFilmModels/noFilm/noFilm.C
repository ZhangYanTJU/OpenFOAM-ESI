/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "noFilm.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noFilm, 0);
addToRunTimeSelectionTable(surfaceFilmModel, noFilm, mesh);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

noFilm::noFilm
(
    const word& modelType,
    const fvMesh& mesh,
    const dimensionedVector& g,
    const word& regionType
)
:
    surfaceFilmModel(),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noFilm::~noFilm()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar noFilm::CourantNumber() const
{
    return 0;
}


tmp<volScalarField::Internal> noFilm::Srho() const
{
    return volScalarField::Internal::New
    (
        IOobject::scopedName("noFilm", "Srho"),
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar(dimMass/dimVolume/dimTime, Zero)
    );
}


tmp<volScalarField::Internal> noFilm::Srho(const label i) const
{
    return volScalarField::Internal::New
    (
        IOobject::scopedName("noFilm", "Srho(" + Foam::name(i) + ")"),
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar(dimMass/dimVolume/dimTime, Zero)
    );
}


tmp<volScalarField::Internal> noFilm::Sh() const
{
    return volScalarField::Internal::New
    (
        IOobject::scopedName("noFilm", "Sh"),
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
    );
}


void noFilm::evolve()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
