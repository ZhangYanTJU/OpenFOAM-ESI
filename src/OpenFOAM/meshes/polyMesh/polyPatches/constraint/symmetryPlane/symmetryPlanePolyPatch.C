/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
    Copyright (C) 2022,2024 OpenCFD Ltd.
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

#include "symmetryPlanePolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "symmetryPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(symmetryPlanePolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, symmetryPlanePolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, symmetryPlanePolyPatch, dictionary);
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

void Foam::symmetryPlanePolyPatch::calcGeometry(PstreamBuffers&)
{
    if (n_ == vector::rootMax)
    {
        if (returnReduceOr(size()))
        {
            // Instead of using the average unit-normal use an area weighted
            // average instead. This avoids problem when adding zero-sized
            // faces since these have a calculated area vector of (0 0 0)
            const auto& areas = faceAreas();

            n_ = gSum(areas).normalise(ROOTVSMALL);

            if (debug)
            {
                Info<< "Patch " << name() << " calculated average normal "
                    << n_ << endl;
            }


            // Check the symmetry plane is planar
            forAll(areas, facei)
            {
                const scalar a = mag(areas[facei]);

                // Calculate only if non-zero area
                if (a > ROOTVSMALL)
                {
                    const vector nf(areas[facei]/a);

                    if (magSqr(n_ - nf) > SMALL)
                    {
                        FatalErrorInFunction
                            << "Symmetry plane '" << name()
                            << "' is not planar." << endl
                            << "At local face at "
                            << primitivePatch::faceCentres()[facei]
                            << " the normal " << nf
                            << " differs from the average normal " << n_
                            << " by " << magSqr(n_ - nf) << endl
                            << "Either split the patch into planar parts"
                            << " or use the " << symmetryPolyPatch::typeName
                            << " patch type"
                            << exit(FatalError);
                    }
                }
            }
        }
    }
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::symmetryPlanePolyPatch::symmetryPlanePolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, size, start, index, bm, patchType),
    n_(vector::rootMax)
{}


Foam::symmetryPlanePolyPatch::symmetryPlanePolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, dict, index, bm, patchType),
    n_(vector::rootMax)
{}


Foam::symmetryPlanePolyPatch::symmetryPlanePolyPatch
(
    const symmetryPlanePolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    n_(pp.n_)
{}


Foam::symmetryPlanePolyPatch::symmetryPlanePolyPatch
(
    const symmetryPlanePolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart),
    n_(pp.n_)
{}


Foam::symmetryPlanePolyPatch::symmetryPlanePolyPatch
(
    const symmetryPlanePolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    polyPatch(pp, bm, index, mapAddressing, newStart),
    n_(pp.n_)
{}


// ************************************************************************* //
