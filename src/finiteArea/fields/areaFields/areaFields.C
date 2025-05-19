/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2018-2025 OpenCFD Ltd.
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

#include "faMesh.H"
#include "areaFields.H"
#include "coupledFaPatchField.H"
#include "registerSwitch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTemplateTypeNameAndDebug(areaScalarField::Internal, 0);
defineTemplateTypeNameAndDebug(areaVectorField::Internal, 0);
defineTemplateTypeNameAndDebug(areaSphericalTensorField::Internal, 0);
defineTemplateTypeNameAndDebug(areaSymmTensorField::Internal, 0);
defineTemplateTypeNameAndDebug(areaTensorField::Internal, 0);

defineTemplateTypeNameAndDebug(areaScalarField, 0);
defineTemplateTypeNameAndDebug(areaVectorField, 0);
defineTemplateTypeNameAndDebug(areaSphericalTensorField, 0);
defineTemplateTypeNameAndDebug(areaSymmTensorField, 0);
defineTemplateTypeNameAndDebug(areaTensorField, 0);

defineTemplateDebugSwitchWithName
(
    areaScalarField::Boundary,
    "areaScalarField::Boundary",
    0
);
defineTemplateDebugSwitchWithName
(
    areaVectorField::Boundary,
    "areaVectorField::Boundary",
    0
);
defineTemplateDebugSwitchWithName
(
    areaSphericalTensorField::Boundary,
    "areaSphericalTensorField::Boundary",
    0
);
defineTemplateDebugSwitchWithName
(
    areaSymmTensorField::Boundary,
    "areaSymmTensorField::Boundary",
    0
);
defineTemplateDebugSwitchWithName
(
    areaTensorField::Boundary,
    "areaTensorField::Boundary",
    0
);

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Specializations

namespace Foam
{

template<>
tmp<GeometricField<scalar, faPatchField, areaMesh>>
GeometricField<scalar, faPatchField, areaMesh>::component
(
    const direction
) const
{
    return *this;
}

template<>
void GeometricField<scalar, faPatchField, areaMesh>::replace
(
    const direction,
    const GeometricField<scalar, faPatchField, areaMesh>& gsf
)
{
    *this == gsf;
}


#undef  fieldChecks
#define fieldChecks(Type)                                                     \
template<>                                                                    \
bool GeometricBoundaryField<Type, faPatchField, areaMesh>::check() const      \
{                                                                             \
    return checkConsistency<coupledFaPatchField<Type>>                        \
    (                                                                         \
        FieldBase::localBoundaryTolerance_,                                   \
       !(debug&4)  /* make into warning if debug&4 */                         \
    );                                                                        \
}

fieldChecks(scalar);
fieldChecks(vector);
fieldChecks(sphericalTensor);
fieldChecks(symmTensor);
fieldChecks(tensor);

#undef fieldChecks

} // End namespace Foam


// * * * * * * * * * * * * * * * * Global Data * * * * * * * * * * * * * * * //

// Note hard-coded values are more reliable than other alternatives

const Foam::wordList Foam::fieldTypes::area
({
    "areaScalarField",
    "areaVectorField",
    "areaSphericalTensorField",
    "areaSymmTensorField",
    "areaTensorField"
});

const Foam::wordList Foam::fieldTypes::area_internal
({
    "areaScalarField::Internal",
    "areaVectorField::Internal",
    "areaSphericalTensorField::Internal",
    "areaSymmTensorField::Internal",
    "areaTensorField::Internal"
});


// ************************************************************************* //
