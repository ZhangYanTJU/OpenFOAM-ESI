/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 PCOpt/NTUA
    Copyright (C) 2022      FOSS GP
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
    compressedGeometricField

Description
    Base class for compression and decompression of the primal fields

SourceFiles
    compressedGeometricField.C

\*---------------------------------------------------------------------------*/


#ifndef compressedGeometricField_H
#define compressedGeometricField_H

#include "GeometricField.H"
#include "fieldTypes.H"
#include "storageParameters.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class compressedGeometricField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class PatchField, class GeoMesh>
class compressedGeometricField
{
    static label counter;

private:

    // Private Member Functions

        //- No copy construct
        compressedGeometricField(const compressedGeometricField&) = delete;

        //- No copy assignment
        void operator=(const compressedGeometricField&) = delete;


protected:

    // Protected data

        //- Reference to the field to be compressed,
        //- i.e. the field used from primalSolver
        GeometricField<Type, PatchField, GeoMesh>& field_;

        //- Name of the field to be compressed
        const label iPtr_;

        const dictionary dict_;
        const fvMesh& mesh_;

        //- Is the first non-empty direction
        label validCompon_;

        //- Current S/N of the time-step that is compressed in the present
        //- object
        const label& k_;

        //- Reference to the compression parameters
        storageParameters& storageParams_;

        //- Are all unnecessary actions are avoided to benchmark the code?
        const bool timing_;

        //- The vector of solved-for directions in mesh.
        //  1 indicates valid direction and -1 an invalid direction.
        const Vector<label>& solDirs_;

        //- Per patch, whether to store/compress fields using the
        //- appropriate algorithm.
        //  Patches with uniform BCs, if a single value is chosen to be
        //  stored for them, have a false value
        boolList toBeStored_;

        //- Per patch of a surfaceScalarField, whether they must be zeroed
        //- upon retrieval.
        //  To do so a zero value is stored for each of them.
        boolList surfacePatchesToBeZeroed_;

        //- Per patch, whether to store/compress fields using the
        //- appropriate algorithm, irrespective of storing a uniform or
        //- non-uniform value
        boolList toBeStoredUniformOrNot_;

        //- Threshold to avoid division with Zero
        scalar eps_;

        //- List holding all of the storageMetrics information
        scalarList storageMetrics_;

        //- Initial memory size of the field
        scalar& initialSize_;

        //- Initial memory size of the part of the field, which will be stored
        //  (e.g. excluding the coupled patches)
        scalar& uncompressedRoughSize_;

        //- Initial memory size of the part of the field that is sent to Zfp
        //- compressor.
        //  It may be larger than uncompressedRoughSize_, if
        //  rearrangementType=2 and meshMode_=1. Needed to compute the pure
        //  compression ratio of zfp.
        scalar& uncompressedSize_;

        //- Compressed memory size of the field
        scalar& compressedSize_;

        //- Total size of elements at all lists to be compressed
        label totalSize_;


    // Protected Member Functions

        //- Fill in uniformBoundaryValuePtr_ with the uniform boundary field
        //- values
        virtual void storeUniformBoundaryValues() = 0;

        //- Restore the uniform boundary field values based on what is stored
        //- in uniformBoundaryValuePtr_
        virtual void restoreUniformBoundaryValues
        (
            GeometricField<Type, PatchField, GeoMesh>& field
        ) = 0;

        //- Const and non-const function to return a reference to
        //- field[pos][componI] for a vectorField and to field[pos] for a
        //- scalarField
        inline scalar& getValue
        (
            Field<Type>& field,
            const label& pos,
            const label& componI
        ) const;

        inline const scalar& getValue
        (
            const Field<Type>& field,
            const label& pos,
            const label& componI
        ) const;

        //- Const and non-const function to return a reference to
        //- field[pos][componI] for a List<vector> and to field[pos] for a
        //- List<scalar>
        inline scalar& getValue
        (
            List<Type>& field,
            const label& pos,
            const label& componI
        ) const;

        inline const scalar& getValue
        (
            const List<Type>& field,
            const label& pos,
            const label& componI
        ) const;

        //- Fill in toBeSorted_ which defines if the scalarField component of
        //- a boundary condition must be stored or not, per patch
        //  (empty, coupled and zero sized patches do not get stored)
        virtual void toBeStored();

        //- Since partial template specialization is not allowed,
        //- define two different functions with the same arguments
        //- and specialize checks there
        void checkFvsPatchFieldStorage
        (
            const PatchField<Type>& patchField,
            const label pI
        );

        void checkFvPatchFieldStorage
        (
            const PatchField<Type>& patchField,
            const label patchI
        );

        //- Compute initial size of the field to be compressed
        void computeSize();

        //- Construct a list containing all of the compression metrics
        void gatherStorageMetrics();

        //- Set a valid solution direction
        void setValidComponent();

        //- Called before correcting the BCs
        //  The function is specialized for the different geometric fields.
        inline void handleParticularities
        (
            GeometricField<Type, PatchField, GeoMesh>& field
        );

        //- Correct boundary conditions of coupled patches.
        //  The function is specialized for the different geometric fields.
        inline void correctBCs
        (
            GeometricField<Type, PatchField, GeoMesh>& field
        );


public:

    //- Runtime type information
    TypeName("compressedGeometricField");

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            compressedGeometricField,
            dictionary,
            (
                GeometricField<Type, PatchField, GeoMesh>& field,
                storageParameters& storageParams,
                const label iPtr,
                const label& k
            ),
            (field, storageParams, iPtr, k)
        );


    // Constructors
    compressedGeometricField
    (
        GeometricField<Type, PatchField, GeoMesh>& field,
        storageParameters& storageParams,
        const label iPtr,
        const label& k
    );


    // Destructors
    virtual ~compressedGeometricField() = default;


    // Selectors

        //- Return a pointer to a new compressedGeometricField
        static autoPtr<compressedGeometricField<Type, PatchField, GeoMesh>> New
        (
            GeometricField<Type, PatchField, GeoMesh>& field,
            storageParameters& storageParams,
            const label iPtr,
            const label& k
        );


    // Member Functions

        //- Compress the field
        virtual void compress();

        //- Decompress the field
        virtual void decompress
        (
            GeometricField<Type, PatchField, GeoMesh>& field
        );

        void decompress();

        const word& name() const;

        //- Get a const reference to the compression metrics
        virtual const scalarList& storageMetrics() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "compressedGeometricField.C"
#endif


// Add the patch constructor functions to the hash tables

#define makeCompressedGeometricTypeField(SS, Type, PatchField, GeoMesh)        \
    typedef SS<Type, PatchField, GeoMesh>                                      \
        SS##Type##PatchField##GeoMesh##_;                                      \
                                                                               \
    defineNamedTemplateTypeNameAndDebug                                        \
        (SS##Type##PatchField##GeoMesh##_, 0);                                 \
                                                                               \
    compressedGeometricField<Type, PatchField, GeoMesh>                        \
    ::adddictionaryConstructorToTable<SS<Type, PatchField, GeoMesh>>           \
        add##SS##Type##PatchField##GeoMesh##dictionaryConstructorToTable_;     \

#define makeCompressedGeometricTypeFieldNoAdd(SS, Type, PatchField, GeoMesh)   \
    typedef SS<Type, PatchField, GeoMesh>                                      \
        SS##Type##PatchField##GeoMesh##_;                                      \
                                                                               \
    defineNamedTemplateTypeNameAndDebug                                        \
        (SS##Type##PatchField##GeoMesh##_, 0);

#define makeCompressedGeometricField(SS)                                       \
                                                                               \
makeCompressedGeometricTypeField(SS, scalar, fvPatchField,  volMesh)           \
makeCompressedGeometricTypeField(SS, scalar, fvsPatchField, surfaceMesh)       \
makeCompressedGeometricTypeField(SS, vector, fvPatchField,  volMesh)

#define makeCompressedGeometricFieldNoAdd(SS)                                  \
                                                                               \
makeCompressedGeometricTypeFieldNoAdd(SS, scalar, fvPatchField,  volMesh)      \
makeCompressedGeometricTypeFieldNoAdd(SS, scalar, fvsPatchField, surfaceMesh)  \
makeCompressedGeometricTypeFieldNoAdd(SS, vector, fvPatchField,  volMesh)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
