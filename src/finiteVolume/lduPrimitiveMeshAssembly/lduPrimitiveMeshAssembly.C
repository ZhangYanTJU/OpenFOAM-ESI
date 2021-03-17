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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(lduPrimitiveMeshAssembly, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::lduPrimitiveMeshAssembly::totalSize
(
    const UPtrList<fvMesh>& meshes
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
    const UPtrList<fvMesh>& meshes
)
{
    if (meshes.size() == 1)
    {
        return 0;
    }
    else
    {
        forAll(meshes, meshi)
        {
            if (meshes[meshi].name() == pp.boundaryMesh().mesh().name())
            {
                return meshi;
            }
        }
    }
    return -1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lduPrimitiveMeshAssembly::lduPrimitiveMeshAssembly
(
    const UPtrList<fvMesh>& meshes,
    const IOobject& io
)
:
    regIOobject(io),
    lduPrimitiveMesh(totalSize(meshes)),
    meshes_(meshes),
    timeIndex_(meshes[0].time().timeIndex())
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
    magSfFaceBoundMap_.setSize(nMeshes);


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
    const fvMesh& mesh,
    const IOobject& io
)
:
    regIOobject(io),
    lduPrimitiveMesh(mesh.lduAddr().size()),
    meshes_(1),
    timeIndex_(mesh.time().timeIndex())
{
    //meshes_.setSize(1);
    meshes_.set(0, const_cast<fvMesh*>(&mesh));
    const label nMeshes(1);
    patchMap_.setSize(nMeshes);
    faceBoundMap_.setSize(nMeshes);
    cellBoundMap_.setSize(nMeshes);
    faceMap_.setSize(nMeshes);
    patchLocalToGlobalMap_.setSize(nMeshes);
    magSfFaceBoundMap_.setSize(nMeshes);

    // Determine cellOffset and faceOffset
    cellOffsets_.setSize(nMeshes);
    cellOffsets_[0] = 0;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::lduPrimitiveMeshAssembly::fluxRequired(const word& name) const
{
    bool flux = false;
    for (label i=0; i < meshes_.size(); ++i)
    {
        if (meshes_[i].fluxRequired(name))
        {
            flux = true;
            break;
        }
    }
    return flux;

}

// ************************************************************************* //
