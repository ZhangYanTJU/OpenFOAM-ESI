/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

#include "processorCyclicPointPatch.H"
#include "pointBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(processorCyclicPointPatch, 0);

addToRunTimeSelectionTable
(
    facePointPatch,
    processorCyclicPointPatch,
    polyPatch
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

processorCyclicPointPatch::processorCyclicPointPatch
(
    const polyPatch& patch,
    const pointBoundaryMesh& bm
)
:
    processorPointPatch(patch, bm),
    procCycPolyPatch_(refCast<const processorCyclicPolyPatch>(patch))
{}


Foam::processorCyclicPointPatch::processorCyclicPointPatch
(
    const processorCyclicPointPatch& patch,
    const pointBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const labelUList& reversePointMap
)
:
    processorPointPatch(patch, bm, index, mapAddressing, reversePointMap),
    procCycPolyPatch_(refCast<const processorCyclicPolyPatch>(patch.patch()))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

processorCyclicPointPatch::~processorCyclicPointPatch()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
