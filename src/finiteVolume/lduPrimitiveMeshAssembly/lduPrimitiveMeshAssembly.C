/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "lduPrimitiveMeshAssembly.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(lduPrimitiveMeshAssembly, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::lduPrimitiveMeshAssembly::totalSize
(
    const UPtrList<lduMesh>& meshes
)
{
    label tot = 0;

    forAll(meshes, meshi)
    {
        tot += meshes[meshi].lduAddr().size();
    }

    return tot;
}


Foam::label Foam::lduPrimitiveMeshAssembly::findNbrMeshId
(
    const polyPatch& pp,
    const label iMesh
) const
{
    if (pp.neighbRegionID() != "none")
    {
        forAll(meshes_, meshi)
        {
            if (meshes_[meshi].thisDb().name() == pp.neighbRegionID())
            {
                return meshi;
            }
        }
    }
    else
    {
        return iMesh;
    }
    return -1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lduPrimitiveMeshAssembly::lduPrimitiveMeshAssembly
(
    const UPtrList<lduMesh>& meshes,
    const IOobject& io
)
:
    regIOobject(io),
    lduPrimitiveMesh(totalSize(meshes)),
    meshes_(meshes),
    timeIndex_(io.time().timeIndex())
{
    forAll(meshes, meshi)
    {
        if (meshes[meshi].comm() != comm())
        {
            WarningInFunction
                << "Communicator " << meshes[meshi].comm()
                << " at index " << meshi
                << " differs between meshes " << nl;
        }
    }

    const label nMeshes = meshes.size();
    patchMap_.setSize(nMeshes);
    faceBoundMap_.setSize(nMeshes);
    cellBoundMap_.setSize(nMeshes);
    faceMap_.setSize(nMeshes);
    patchLocalToGlobalMap_.setSize(nMeshes);
    facePatchFaceMap_.setSize(nMeshes);

    // Determine cellOffset and faceOffset
    cellOffsets_.setSize(1+nMeshes);
    cellOffsets_[0] = 0;
    for (label meshi=0; meshi < nMeshes; ++meshi)
    {
        cellOffsets_[meshi+1] =
            cellOffsets_[meshi] + meshes[meshi].lduAddr().size();
    }
}


Foam::lduPrimitiveMeshAssembly::lduPrimitiveMeshAssembly
(
    const lduMesh& mesh,
    const IOobject& io
)
:
    regIOobject(io),
    lduPrimitiveMesh(mesh.lduAddr().size()),
    meshes_(1),
    timeIndex_(io.time().timeIndex())
{
    meshes_.set(0, const_cast<lduMesh*>(&mesh));
    const label nMeshes(1);
    patchMap_.setSize(nMeshes);
    faceBoundMap_.setSize(nMeshes);
    cellBoundMap_.setSize(nMeshes);
    faceMap_.setSize(nMeshes);
    patchLocalToGlobalMap_.setSize(nMeshes);
    facePatchFaceMap_.setSize(nMeshes);

    // Determine cellOffset and faceOffset
    cellOffsets_.setSize(nMeshes);
    cellOffsets_[0] = 0;
}


// ************************************************************************* //
