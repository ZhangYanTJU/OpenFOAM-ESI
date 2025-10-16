/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenCFD Ltd.
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

#include "cornerDetectionModel.H"
#include "faMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cornerDetectionModel, 0);
    defineRunTimeSelectionTable(cornerDetectionModel, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::cornerDetectionModel::dihedralAngle
(
    const vector& n0,
    const vector& n1
) const
{
    if (mag(n0) <= VSMALL || mag(n1) <= VSMALL)
    {
#ifdef FULL_DEBUG
        WarningInFunction
            << "Degenerate face normal magnitude (|n| ~ 0). "
            << "Returning 0 for dihedral angle." << nl;
#endif
        return 0;
    }

    const scalar a = mag(n1 - n0);
    const scalar b = mag(n1 + n0);

    // The dihedral angle is calculated as 2*atan2(|n1 - n0|, |n1 + n0|),
    // which gives the angle between the two normals n0 and n1.
    scalar phi = scalar(2)*std::atan2(a, b);

    // Clamp to [0, pi]
    phi = max(0, min(constant::mathematical::pi, phi));

    if (!std::isfinite(phi))
    {
#ifdef FULL_DEBUG
        WarningInFunction
            << "Non-finite dihedral angle computed. Returning 0." << nl;
#endif
        return 0;
    }

    return phi;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cornerDetectionModel::cornerDetectionModel
(
    const faMesh& mesh,
    const regionModels::areaSurfaceFilmModels::liquidFilmBase& film,
    const dictionary& dict
)
:
    mesh_(mesh),
    film_(film),
    dict_(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cornerDetectionModel::~cornerDetectionModel()
{}  // faMesh was forward declared


// ************************************************************************* //
