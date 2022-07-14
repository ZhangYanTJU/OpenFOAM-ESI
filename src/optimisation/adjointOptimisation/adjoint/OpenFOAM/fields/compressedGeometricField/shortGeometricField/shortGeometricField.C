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

#include "shortGeometricField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void
shortGeometricField<Type, PatchField, GeoMesh>::storeUniformBoundaryValues()
{
    forAll(uniformBoundaryValuePtr_, pI)
    {
        const Field<Type>& boundaryField = this->field_.boundaryField()[pI];
        if (this->toBeStoredUniformOrNot_[pI] && !this->toBeStored_[pI])
        {
            if
            (
                std::is_same<PatchField<Type>, fvsPatchField<Type>>::value
            &&  this->surfacePatchesToBeZeroed_[pI]
            )
            {
                // At this point, the maximum magnitude of the boundaryField is
                // for sure below the prescribed threshold
                uniformBoundaryValuePtr_.set(pI, new Type(Zero));
            }
            else
            {
                uniformBoundaryValuePtr_.set
                (
                    pI, new Type(boundaryField.first())
                );
            }
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void
shortGeometricField<Type, PatchField, GeoMesh>::restoreUniformBoundaryValues
(
    GeometricField<Type, PatchField, GeoMesh>& field
)
{
    forAll(uniformBoundaryValuePtr_, pI)
    {
        Field<Type>& boundaryField = field.boundaryFieldRef()[pI];
        if (uniformBoundaryValuePtr_(pI))
        {
            boundaryField = uniformBoundaryValuePtr_[pI];
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void shortGeometricField<Type, PatchField, GeoMesh>::initialize()
{
    this->toBeStored();
    this->computeSize();
    this->uncompressedSize_ = this->uncompressedRoughSize_;
    this->compressedSize_ = this->uncompressedRoughSize_;
    setBuffers();
    storeUniformBoundaryValues();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void shortGeometricField<Type, PatchField, GeoMesh>::setBuffers()
{
    const Field<Type>& internalField = this->field_.primitiveField();
    forAll(internalBufferPtr_, cmptI)
    {
        if (this->solDirs_[cmptI])
        {
            internalBufferPtr_[cmptI].setSize(internalField.size());
        }
    }
    forAll (boundaryBufferPtr_, pI)
    {
        if (this->toBeStored_[pI])
        {
            boundaryBufferPtr_[pI].setSize(pTraits<Type>::nComponents);
            const Field<Type>& boundaryField = this->field_.boundaryField()[pI];
            forAll(boundaryBufferPtr_[pI], cmptI)
            {
                if (this->solDirs_[cmptI])
                {
                    boundaryBufferPtr_[pI][cmptI].setSize(boundaryField.size());
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
shortGeometricField<Type, PatchField, GeoMesh>::shortGeometricField
(
    GeometricField<Type, PatchField, GeoMesh>& field,
    storageParameters& storageParams,
    const label iPtr,
    const label& k
):
    compressedGeometricField<Type, PatchField, GeoMesh>
    (
        field,
        storageParams,
        iPtr,
        k
    ),
    internalBufferPtr_(pTraits<Type>::nComponents),
    boundaryBufferPtr_(field.mesh().boundary().size()),
    uniformBoundaryValuePtr_(field.mesh().boundary().size())
{
    initialize();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void shortGeometricField<Type, PatchField, GeoMesh>::compress()
{
    const Field<Type>& internalField = this->field_.primitiveField();
    forAll(internalBufferPtr_, cmptI)
    {
        if (this->solDirs_[cmptI])
        {
            forAll(internalBufferPtr_[cmptI], iPtr)
            {
                internalBufferPtr_[cmptI][iPtr] =
                    this->getValue(internalField, iPtr, cmptI);
            }
        }
    }
    forAll(boundaryBufferPtr_, pI)
    {
        if (this->toBeStored_[pI])
        {
            const Field<Type>& boundaryField = this->field_.boundaryField()[pI];
            for (label cmptI = 0; cmptI < pTraits<Type>::nComponents; ++cmptI)
            {
                if (this->solDirs_[cmptI])
                {
                    forAll(boundaryBufferPtr_[pI][cmptI], iPtr)
                    {
                        boundaryBufferPtr_[pI][cmptI][iPtr] =
                            this->getValue(boundaryField, iPtr, cmptI);
                    }
                }
            }
        }
    }
    this->gatherStorageMetrics();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void shortGeometricField<Type, PatchField, GeoMesh>::decompress
(
    GeometricField<Type, PatchField, GeoMesh>& field
)
{
    Field<Type>& internalField = field.primitiveFieldRef();
    forAll(internalBufferPtr_, cmptI)
    {
        if (this->solDirs_[cmptI])
        {
            forAll(internalBufferPtr_[cmptI], iPtr)
            {
                scalar& val = this->getValue(internalField, iPtr, cmptI);
                val = internalBufferPtr_[cmptI][iPtr];
            }
        }
        /*else
        {
            forAll(internalBufferPtr_[cmptI], iPtr)
            {
                scalar& val = this->getValue(internalField, iPtr, cmptI);
                val = 0;
            }
        }*/
    }
    forAll(boundaryBufferPtr_, pI)
    {
        if (this->toBeStored_[pI])
        {
            Field<Type>& boundaryField = field.boundaryFieldRef()[pI];
            for (label cmptI = 0; cmptI < pTraits<Type>::nComponents; ++cmptI)
            {
                if (this->solDirs_[cmptI])
                {
                    forAll(boundaryBufferPtr_[pI][cmptI], iPtr)
                    {
                        scalar& val =
                            this->getValue(boundaryField, iPtr, cmptI);
                        val = boundaryBufferPtr_[pI][cmptI][iPtr];
                    }
                }
                /*else
                {
                    forAll (boundaryBufferPtr_[pI], iPtr)
                    {
                        scalar& val = this->getValue(boundaryField, iPtr, cmptI);
                        val = 0;
                    }
                }*/
            }
        }
    }
    this->restoreUniformBoundaryValues(field);
    this->correctBCs(field);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
