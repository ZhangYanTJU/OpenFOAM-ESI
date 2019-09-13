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

Class
    Foam::lduPrimitiveMeshAssembly

Description
    An assembly of lduMatrix that is specific inter-region coupling
    through mapped patches.

SourceFiles
    lduPrimitiveMeshAssembly.C

\*---------------------------------------------------------------------------*/

#ifndef lduPrimitiveMeshAssembly_H
#define lduPrimitiveMeshAssembly_H

#include "fvMesh.H"
#include "lduPrimitiveMesh.H"
#include "mappedPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class lduPrimitiveMeshAssembly Declaration
\*---------------------------------------------------------------------------*/

class lduPrimitiveMeshAssembly
:
    public objectRegistry,
    public lduPrimitiveMesh
{
    // Private Data

        List<labelList> patchAddr_;

        //- Map patches from local to global
        labelListList patchMap_;

        //- Map faces from local to global
        labelListList faceMap_;

        //- Map path face from local to global
        labelListListList faceBoundMap_;

        //- Offset cells
        labelList cellOffsets_;


    // Private Member Functions

        //- Get size of all meshes
        static label totalSize(const UPtrList<fvMesh>&);

        //- Find nrb mesh Id for mapped patches
        static label findNbrMeshId
        (
            const mappedPatchBase& pp,
            const UPtrList<fvMesh>& meshes
        );


        //- No copy construct
        lduPrimitiveMeshAssembly(const lduPrimitiveMeshAssembly&) = delete;

        //- No copy assignment
        void operator=(const lduPrimitiveMeshAssembly&) = delete;


public:

    // Static Data

        // Declare name of the class and its debug switch
        ClassName("lduPrimitiveMeshAssembly");


    // Constructors

        //- Construct from fvMeshes combining mapped faces into internal faces
        lduPrimitiveMeshAssembly
        (
            const UPtrList<fvMesh>& meshes,
            const IOobject& io
        );


    //- Destructor
    virtual ~lduPrimitiveMeshAssembly() = default;


    // Member Functions

    // Access

        //- Return patch addressing
        virtual const labelUList& patchAddr(const label patchNo) const
        {
            return patchAddr_[patchNo];
        }

        //- Return the object registry
        virtual const objectRegistry& thisDb() const
        {
            return *this;
        }

        //- Return patchMap
        const labelListList& patchMap() const
        {
            return patchMap_;
        }

        //- Return boundary face map
        const labelListListList& faceBoundMap() const
        {
            return faceBoundMap_;
        }

        //- Return cellOffsets [1+nMeshes]
        const labelList& cellOffsets() const
        {
            return cellOffsets_;
        }

        //- Return faceMap
        const labelListList& faceMap() const
        {
            return faceMap_;
        }

        //- Return patchAddr
        const List<labelList>& patchAddr() const
        {
            return patchAddr_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
