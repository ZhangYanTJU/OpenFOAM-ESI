/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class UnaryPredicate>
Foam::label Foam::polyBoundaryMesh::nFaces_if(UnaryPredicate pred) const
{
    const polyBoundaryMesh& patches = *this;

    label count = 0;

    for (const polyPatch& pp : patches)
    {
        if (pred(pp))
        {
            count += pp.size();
        }
    }

    return count;
}


template<class UnaryPredicate>
Foam::labelList Foam::polyBoundaryMesh::indices_if(UnaryPredicate pred) const
{
    const polyBoundaryMesh& patches = *this;
    const label total = patches.size();

    labelList patchIDs(total);

    label count = 0;

    for (label patchi = 0; patchi < total; ++patchi)
    {
        if (pred(patches[patchi]))
        {
            patchIDs[count] = patchi;
            ++count;
        }
    }

    patchIDs.resize(count);

    return patchIDs;
}


template<class PatchType>
Foam::labelList Foam::polyBoundaryMesh::indices_if() const
{
    const polyBoundaryMesh& patches = *this;
    const label total = patches.size();

    labelList patchIDs(total);

    label count = 0;

    for (label patchi = 0; patchi < total; ++patchi)
    {
        if (isA<PatchType>(patches[patchi]))
        {
            patchIDs[count] = patchi;
            ++count;
        }
    }

    patchIDs.resize(count);

    return patchIDs;
}


template<class PatchType>
Foam::labelHashSet Foam::polyBoundaryMesh::findPatchIDs() const
{
    const polyBoundaryMesh& patches = *this;
    const label total = patches.size();

    labelHashSet patchIDs;
    patchIDs.reserve(total);

    for (label patchi = 0; patchi < total; ++patchi)
    {
        if (isA<PatchType>(patches[patchi]))
        {
            patchIDs.insert(patchi);
        }
    }

    return patchIDs;
}


// ************************************************************************* //
