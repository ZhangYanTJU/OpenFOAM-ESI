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

\*---------------------------------------------------------------------------*/

#include "wallFvPatch.H"
#include "symmetryPlaneFvPatch.H"
#include "emptyFvsPatchField.H"
#include "processorFvsPatchField.H"
#include "coupledFvPatchField.H"
#include "emptyFvPatchField.H"
#include "zeroGradientFvPatchField.H"
#include "inletOutletFvPatchField.H"
#include "outletInletFvPatchField.H"
#include "freestreamFvPatchField.H"
#include "slipFvPatchField.H"
#include "symmetryFvPatchField.H"
#include "symmetryPlaneFvPatchField.H"
#include "cyclicACMIFvPatchField.H"
#include "bound.H"
#include "compressedGeometricField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
label compressedGeometricField<Type, PatchField, GeoMesh>::counter = 0;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<>
inline scalar&
compressedGeometricField<scalar, fvPatchField, volMesh>::getValue
(
    Field<scalar>& field,
    const label& pos,
    const label& componI
) const
{
    return field[pos];
}


template<>
inline scalar&
compressedGeometricField<scalar, fvsPatchField, surfaceMesh>::getValue
(
    Field<scalar>& field,
    const label& pos,
    const label& componI
) const
{
    return field[pos];
}


template<>
inline scalar&
compressedGeometricField<vector, fvPatchField, volMesh>::getValue
(
    Field<vector>& field,
    const label& pos,
    const label& componI
) const
{
    return field[pos][componI];
}


template<>
inline const scalar&
compressedGeometricField<scalar, fvPatchField, volMesh>::getValue
(
    const Field<scalar>& field,
    const label& pos,
    const label& componI
) const
{
    return field[pos];
}


template<>
inline const scalar&
compressedGeometricField<scalar, fvsPatchField, surfaceMesh>::getValue
(
    const Field<scalar>& field,
    const label& pos,
    const label& componI
) const
{
    return field[pos];
}


template<>
inline const scalar&
compressedGeometricField<vector, fvPatchField, volMesh>::getValue
(
    const Field<vector>& field,
    const label& pos,
    const label& componI
) const
{
    return field[pos][componI];
}


template<>
inline scalar&
compressedGeometricField<scalar, fvPatchField, volMesh>::getValue
(
    List<scalar>& field,
    const label& pos,
    const label& componI
) const
{
    return field[pos];
}


template<>
inline scalar&
compressedGeometricField<scalar, fvsPatchField, surfaceMesh>::getValue
(
    List<scalar>& field,
    const label& pos,
    const label& componI
) const
{
    return field[pos];
}


template<>
inline scalar&
compressedGeometricField<vector, fvPatchField, volMesh>::getValue
(
    List<vector>& field,
    const label& pos,
    const label& componI
) const
{
    return field[pos][componI];
}


template<>
inline const
scalar& compressedGeometricField<scalar, fvPatchField, volMesh>::getValue
(
    const List<scalar>& field,
    const label& pos,
    const label& componI
) const
{
    return field[pos];
}


template<>
inline const
scalar& compressedGeometricField<scalar, fvsPatchField, surfaceMesh>::getValue
(
    const List<scalar>& field,
    const label& pos,
    const label& componI
) const
{
    return field[pos];
}


template<>
inline const
scalar& compressedGeometricField<vector, fvPatchField, volMesh>::getValue
(
    const List<vector>& field,
    const label& pos,
    const label& componI
) const
{
    return field[pos][componI];
}


template<class Type, template<class> class PatchField, class GeoMesh>
void compressedGeometricField<Type, PatchField, GeoMesh>::toBeStored()
{
    setValidComponent();
    toBeStored_.setSize(mesh_.boundary().size(), false);
    toBeStoredUniformOrNot_.setSize(mesh_.boundary().size(), false);
    surfacePatchesToBeZeroed_.setSize(mesh_.boundary().size(), false);
    typedef typename GeometricField<Type, PatchField, GeoMesh>::Boundary
        Boundary;
    const Boundary& boundaryField = field_.boundaryField();
    forAll(mesh_.boundary(), pI)
    {
        const PatchField<Type>& pf = boundaryField[pI];
        if (storageParams_.storeAllBoundaries())
        {
            if (pf.size() != 0)
            {
                // Warning: Possible error upon compression of boundaries
                // with very small values when single precision is used.
                // E.g. flux phi on boundaries with rotatingWallVelocity
                // for the velocity
                toBeStoredUniformOrNot_[pI] = true;
                toBeStored_[pI] = true;
            }
        }
        else
        {
            checkFvsPatchFieldStorage(pf, pI);
            checkFvPatchFieldStorage(pf, pI);
        }
        if (debug > 1)
        {
            Info<< "Field '" << field_.name() << "': "
                << "Storing boundary patch " << pI
                << " '" << mesh_.boundary()[pI].name() << "' "
                << "of type '" << boundaryField.types()[pI] << "'? "
                << Switch(toBeStored_[pI]);
            if (toBeStoredUniformOrNot_[pI] && !toBeStored_[pI])
            {
                Info << "  --> Uniform boundary. A single value is stored.";
            }
            Info << endl;
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void
compressedGeometricField<Type, PatchField, GeoMesh>::checkFvsPatchFieldStorage
(
    const PatchField<Type>& patchField,
    const label pI
)
{
    if
    (
        isA<fvsPatchField<Type>>(patchField)
     && !isA<emptyFvsPatchField<Type>>(patchField)
     && patchField.size()
    )
    {
        toBeStoredUniformOrNot_[pI] = true;
        // Store processor boundaries even if they are not uniform.
        // Not sure what will happen at a latter time-step for
        // compression techniques using time-windows, such as iPGD
        // Processor check only valid for surfaceFields since volFields
        // wont enter here due to !isA<coupledFvPatch>
        if
        (
            !patchField.uniform()
         || storageParams_.storeUniformBoundaries()
         || isA<processorFvsPatchField<Type>>(patchField)
        )
        {
            toBeStored_[pI] = true;
        }
        //  At boundaries with a near-zero value of flux, store a
        //  (single) zero value for phi. This was added after ZFP
        //  crashed trying to compress single-precision values of the
        //  iPGD modes with magnitude below 5.e-30. That problem does
        //  not arise if these values are represented in
        //  double-precision. In case that the maximum flux magnitude
        //  is lower than 1.e-15, it is assumed that the field should
        //  be zero and the non-zero values arise due to arithmetic
        //  error. This threshold has been arbitrarily chosen.
        //  Target  patch types (although not checked in the code):
        //  wall, symmetryPlane

        //  surfacePatchesToBeZeroed becomes true only in case that phi at the
        //  corresponding patch must be stored (toBeStoredUniformOrNot = true)
        //  and is not uniform (patchField.uniform() = false)
        if (toBeStored_[pI])
        {
            // Do the magnitude check here, to avoid unnecessary max
            // mag calculations
            if (max(mag(patchField)) < 1.e-15)
            {
                toBeStored_[pI] = false;
                surfacePatchesToBeZeroed_[pI] = true;
            }
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void
compressedGeometricField<Type, PatchField, GeoMesh>::checkFvPatchFieldStorage
(
    const PatchField<Type>& patchField,
    const label pI
)
{
    const fvPatch& patch = mesh_.boundary()[pI];
    if
    (
         isA<fvPatchField<Type>>(patchField)
     && !isA<emptyFvPatchField<Type>>(patchField)
     && !isA<coupledFvPatchField<Type>>(patchField)
     && !isA<wallFvPatch>(patch)
     && !isA<cyclicAMIFvPatchField<Type>>(patchField)
     && !isA<zeroGradientFvPatchField<Type>>(patchField)
     && !isA<inletOutletFvPatchField<Type>>(patchField)
     && !isA<outletInletFvPatchField<Type>>(patchField)
     && !isA<freestreamFvPatchField<Type>>(patchField)
     && !isA<slipFvPatchField<Type>>(patchField)
     && !isA<symmetryFvPatchField<Type>>(patchField)
     && !isA<symmetryPlaneFvPatchField<Type>>(patchField)
     &&  patchField.size()
    )
    {
        toBeStoredUniformOrNot_[pI] = true;
        if (!patchField.uniform() || storageParams_.storeUniformBoundaries())
        {
            toBeStored_[pI] = true;
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void compressedGeometricField<Type, PatchField, GeoMesh>::setValidComponent()
{
    validCompon_ = 0;
    if (!solDirs_[0])
    {
        validCompon_ = 1;
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void compressedGeometricField<Type, PatchField, GeoMesh>::computeSize()
{
    initialSize_ = Zero;
    uncompressedRoughSize_ = Zero;
    totalSize_  = Zero;
    const Field<Type>& internalField = field_.primitiveField();
    for (label cmptI = 0; cmptI < pTraits<Type>::nComponents; ++cmptI)
    {
        initialSize_ += internalField.size()*sizeof(scalar);
        if (solDirs_[cmptI])
        {
            uncompressedRoughSize_ += internalField.size()*sizeof(scalar);
            totalSize_ += internalField.size();
        }
        forAll(mesh_.boundary(), pI)
        {
            const Field<Type>& boundaryField = field_.boundaryField()[pI];
            initialSize_ += boundaryField.size()*sizeof(scalar);
            if (toBeStored_[pI] && solDirs_[cmptI])
            {
                uncompressedRoughSize_ += boundaryField.size()*sizeof(scalar);
                totalSize_ += boundaryField.size();
            }
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void compressedGeometricField<Type, PatchField, GeoMesh>::gatherStorageMetrics()
{
    if (!timing_)
    {
        for (scalar& sm : storageMetrics_)
        {
            reduce(sm, sumOp<scalar>());
        }
    }
}


template<>
inline void
compressedGeometricField<scalar, fvPatchField, volMesh>::handleParticularities
(
    volScalarField& field
)
{
    // Due to compression error, turbulence variables might take on
    // negative values. Correct here
    // TODO: A little dangerous way to identify turbuelence fields
    if (iPtr_ == 3 || iPtr_ == 4)
    {
        bound(field, dimensionedScalar(field.dimensions(), Zero));
    }
}


template<>
inline void
compressedGeometricField<vector, fvPatchField, volMesh>::handleParticularities
(
    volVectorField& field
)
{
    vectorField& internalField = field_.primitiveFieldRef();
    const labelList& cellIDs = storageParams_.wallList();
    label cellCounter(0);
    forAll(cellIDs, i)
    {
        vector& vec = internalField[cellIDs[i]];
        if (mag(vec) < SMALL)
        {
            cellCounter++;
            forAll(vec, j)
            {
                if (solDirs_[j])
                {
                    vec[j] = SMALL;
                }
            }
        }
    }
    DebugInfo
        << "Modified " << field.name() << " in " << cellCounter
        << " cell(s) next to wall boundaries" << endl;
}


template<>
inline void
compressedGeometricField<scalar, fvsPatchField, surfaceMesh>::handleParticularities
(
    surfaceScalarField& field
)
{
    // Does nothing
}


template<>
inline void compressedGeometricField<scalar, fvPatchField, volMesh>::correctBCs
(
    volScalarField& field
)
{
    handleParticularities(field);
    if (!storageParams_.storeAllBoundaries())
    {
        field.correctBoundaryConditions();
    }
}


template<>
inline void compressedGeometricField<vector, fvPatchField, volMesh>::correctBCs
(
    volVectorField& field
)
{
    handleParticularities(field);
    if (!storageParams_.storeAllBoundaries())
    {
        field.correctBoundaryConditions();
    }
}


template<>
inline void
compressedGeometricField<scalar, fvsPatchField, surfaceMesh>::correctBCs
(
    surfaceScalarField& field
)
{
    // Does nothing
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
compressedGeometricField<Type, PatchField, GeoMesh>::compressedGeometricField
(
    GeometricField<Type, PatchField, GeoMesh>& field,
    storageParameters& storageParams,
    const label iPtr,
    const label& k
):
    field_(field),
    iPtr_(iPtr),

    dict_(storageParams.dict()),
    mesh_(field.mesh()),
    validCompon_(0),
    k_(k),

    storageParams_(storageParams),
    timing_(storageParams.timing()),
    solDirs_(storageParams.solDirs()),
    eps_(1.e-12),

    storageMetrics_(4, Zero),
    initialSize_(storageMetrics_[0]),
    uncompressedRoughSize_(storageMetrics_[1]),
    uncompressedSize_(storageMetrics_[2]),
    compressedSize_(storageMetrics_[3]),
    totalSize_(Zero)
{
    if (Pstream::master()) counter++;
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
autoPtr<compressedGeometricField<Type, PatchField, GeoMesh>>
compressedGeometricField<Type, PatchField, GeoMesh>::New
(
    GeometricField<Type, PatchField, GeoMesh>& field,
    storageParameters& storageParams,
    const label iPtr,
    const label& k
)
{
    const word& type = storageParams.compressionMethod()[iPtr];
    auto cstrIter = dictionaryConstructorTablePtr_->cfind(type);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(storageParams.dict())
            << "Unknown compressedGeometricField type " << type << nl << nl
            << "Valid compressedGeometricField types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<compressedGeometricField<Type, PatchField, GeoMesh>>
        (cstrIter()( field, storageParams, iPtr, k ));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void compressedGeometricField<Type, PatchField, GeoMesh>::compress()
{
    NotImplemented
}


template<class Type, template<class> class PatchField, class GeoMesh>
void compressedGeometricField<Type, PatchField, GeoMesh>::decompress
(
    GeometricField<Type, PatchField, GeoMesh>& field
)
{
    NotImplemented
}


template<class Type, template<class> class PatchField, class GeoMesh>
void compressedGeometricField<Type, PatchField, GeoMesh>::decompress()
{
    this->decompress(field_);
}


template<class Type, template<class> class PatchField, class GeoMesh>
const word& compressedGeometricField<Type, PatchField, GeoMesh>::name() const
{
    return field_.name();
}


template<class Type, template<class> class PatchField, class GeoMesh>
const scalarList&
compressedGeometricField<Type, PatchField, GeoMesh>::storageMetrics() const
{
    return storageMetrics_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
