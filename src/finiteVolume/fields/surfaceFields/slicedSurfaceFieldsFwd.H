/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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

InClass
    Foam::slicedSurfaceFields

\*---------------------------------------------------------------------------*/

#ifndef Foam_slicedSurfaceFieldsFwd_H
#define Foam_slicedSurfaceFieldsFwd_H

#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template
<
    class Type,
    template<class> class PatchField,
    template<class> class SlicedPatchField,
    class GeoMesh
>
class SlicedGeometricField;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class surfaceMesh;
template<class Type> class fvsPatchField;
template<class Type> class slicedFvsPatchField;


//- A (volume) surface field slice for a given type
template<class Type>
using SurfaceSliceField =
SlicedGeometricField<Type, fvsPatchField, slicedFvsPatchField, surfaceMesh>;


// Typedefs

typedef
    SlicedGeometricField
    <scalar, fvsPatchField, slicedFvsPatchField, surfaceMesh>
    slicedSurfaceScalarField;

typedef
    SlicedGeometricField
    <vector, fvsPatchField, slicedFvsPatchField, surfaceMesh>
    slicedSurfaceVectorField;

typedef
    SlicedGeometricField
    <sphericalTensor, fvsPatchField, slicedFvsPatchField, surfaceMesh>
    slicedSurfaceSphericalTensorField;

typedef
    SlicedGeometricField
    <symmTensor, fvsPatchField, slicedFvsPatchField, surfaceMesh>
    slicedSurfaceSymmTensorField;

typedef
    SlicedGeometricField
    <tensor, fvsPatchField, slicedFvsPatchField, surfaceMesh>
    slicedSurfaceTensorField;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
