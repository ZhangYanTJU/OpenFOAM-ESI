/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021-2024 OpenCFD Ltd.
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

#include "pointFieldDecomposer.H"
#include "processorPointPatchFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh>>
Foam::pointFieldDecomposer::decomposeField
(
    const GeometricField<Type, pointPatchField, pointMesh>& field
) const
{
    // Create a list of pointers for the patchFields
    PtrList<pointPatchField<Type>> patchFields(boundaryAddressing_.size());

    // Create and map the patch field values
    forAll(boundaryAddressing_, patchi)
    {
        if (patchFieldDecomposerPtrs_.set(patchi))
        {
            patchFields.set
            (
                patchi,
                pointPatchField<Type>::New
                (
                    field.boundaryField()[boundaryAddressing_[patchi]],
                    procMesh_.boundary()[patchi],
                    pointPatchField<Type>::Internal::null(),
                    patchFieldDecomposerPtrs_[patchi]
                )
            );
        }
        else
        {
            patchFields.set
            (
                patchi,
                new processorPointPatchField<Type>
                (
                    procMesh_.boundary()[patchi],
                    pointPatchField<Type>::Internal::null()
                )
            );
        }
    }

    // Create the field for the processor
    return GeometricField<Type, pointPatchField, pointMesh>::New
    (
        field.name(),
        IOobject::NO_REGISTER,
        procMesh_,
        field.dimensions(),
        // Internal field - mapped values
        Field<Type>(field.primitiveField(), pointAddressing_),
        // Boundary field
        patchFields
    );
}


template<class GeoField>
void Foam::pointFieldDecomposer::decomposeFields
(
    const PtrList<GeoField>& fields
) const
{
    for (const auto& fld : fields)
    {
        decomposeField(fld)().write();
    }
}


// ************************************************************************* //
