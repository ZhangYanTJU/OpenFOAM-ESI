/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2024 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::rawTopoChangerFvMesh::setUnmappedValues
(
    GeometricField<Type, PatchField, GeoMesh>& fld,
    const bitSet& mappedFace,
    const GeometricField<Type, PatchField, GeoMesh>& baseFld
)
{
    //Pout<< "Checking field " << fld.name() << endl;

    forAll(fld.boundaryField(), patchi)
    {
        auto& fvp = const_cast<PatchField<Type>&>(fld.boundaryField()[patchi]);

        const label start = fvp.patch().start();
        forAll(fvp, i)
        {
            if (!mappedFace[start+i])
            {
                //Pout<< "** Resetting unassigned value on patch "
                //    << fvp.patch().name()
                //    << " localface:" << i
                //    << " to:" << baseFld.boundaryField()[patchi][i] << endl;
                fvp[i] = baseFld.boundaryField()[patchi][i];
            }
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::rawTopoChangerFvMesh::zeroUnmappedValues
(
    const bitSet& mappedFace
) const
{
    typedef GeometricField<Type, PatchField, GeoMesh> FieldType;

    std::unique_ptr<FieldType> zeroFieldPtr;

    for (const word& fldName : names<FieldType>())
    {
        FieldType& fld = lookupObjectRef<FieldType>(fldName);
        //Pout<< "Checking field " << fld.name() << endl;

        if (!zeroFieldPtr)
        {
            zeroFieldPtr = std::make_unique<FieldType>
            (
                this->newIOobject("zero"),
                *this,
                Foam::zero{},
                dimless
            );
        }

        zeroFieldPtr->dimensions().reset(fld.dimensions());

        setUnmappedValues(fld, mappedFace, *zeroFieldPtr);
    }
}


// ************************************************************************* //
