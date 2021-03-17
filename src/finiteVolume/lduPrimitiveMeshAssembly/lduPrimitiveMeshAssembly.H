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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class lduPrimitiveMeshAssembly Declaration
\*---------------------------------------------------------------------------*/
class lduPrimitiveMeshAssembly
:
    public regIOobject,
    public lduPrimitiveMesh
{
public:

private:

    // Private Data

        //- UPtrList for fvMesh
        UPtrList<fvMesh> meshes_;

        //- Overwrite patchAddr from lduPrimitiveMesh for assembly ldu address
        List<labelList> patchAddr_;

        //- Map patches from local (-1 for non-existing)
        labelListList patchMap_;

        //- Map patches from local to global
        labelListList patchLocalToGlobalMap_;

        //- Map internal faces from local to global
        labelListList faceMap_;

        //- Map patch local sub-face to global new internal face
        labelListListList faceBoundMap_;

        //- Map patch local sub-face to nbrCellId
        labelListListList cellBoundMap_;

        //- Map patch local sub-face to local patch face
        labelListListList magSfFaceBoundMap_;

        //- Offset cells
        labelList cellOffsets_;

        //- Time index
        label timeIndex_;


    // Private Member Functions

        //- Get size of all meshes
        static label totalSize(const UPtrList<fvMesh>&);

        //- Find nrb mesh Id for mapped patches
        static label findNbrMeshId
        (
            const polyPatch& pp,
            const UPtrList<fvMesh>& meshes
        );

        //- No copy construct
        lduPrimitiveMeshAssembly(const lduPrimitiveMeshAssembly&) = delete;

        //- No copy assignment
        void operator=(const lduPrimitiveMeshAssembly&) = delete;


public:

    // Static Data

        // Declare name of the class and its debug switch
        TypeName("lduPrimitiveMeshAssembly");


    // Constructors

        //- Construct from fvMeshes combining mapped faces into internal faces
        lduPrimitiveMeshAssembly
        (
            const UPtrList<fvMesh>& meshes,
            const IOobject& io
        );

        //- Construct from one mesh
        lduPrimitiveMeshAssembly
        (
            const fvMesh&,
            const IOobject& io
        );


    // Public functions

        //- Update mappings
        template<class Type>
        void update
        (
            UPtrList<GeometricField<Type, fvPatchField, volMesh>> psis
        );


    //- Destructor
    virtual ~lduPrimitiveMeshAssembly() = default;


    // Member Functions

        //- Time index
        label timeIndex()
        {
            return timeIndex_;
        }

        label timeIndex() const
        {
            return timeIndex_;
        }


        //- Return the object registry
        virtual const objectRegistry& thisDb() const
        {
            return meshes_[0].thisDb();
        }

    // Access

        //- Return patch addressing
        virtual const labelUList& patchAddr(const label patchNo) const
        {
            return patchAddr_[patchNo];
        }


        //- Return patchMap
        const labelListList& patchMap() const
        {
            return patchMap_;
        }

        //- Return patchLocalToGlobalMap
        const labelListList& patchLocalToGlobalMap() const
        {
            return patchLocalToGlobalMap_;
        }

        //- Return boundary face map
        const labelListListList& faceBoundMap() const
        {
            return faceBoundMap_;
        }

        //- Return boundary face cell
        const labelListListList& cellBoundMap() const
        {
            return cellBoundMap_;
        }

        const labelListListList& magSfFaceBoundMap() const
        {
            return magSfFaceBoundMap_;
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

        //- Return fvMeshes
        const UPtrList<fvMesh>& meshes() const
        {
            return meshes_;
        }

        //- Flux required
        bool fluxRequired(const word& name) const;

        //- Write Data
        virtual bool writeData(Ostream&) const
        {
            NotImplemented;
            return false;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "lduPrimitiveMeshAssemblyTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
