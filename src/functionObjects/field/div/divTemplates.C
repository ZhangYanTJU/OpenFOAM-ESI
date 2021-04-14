/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "fvcDiv.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class FieldType>
bool Foam::functionObjects::div::calcDiv()
{
    if (foundObject<FieldType>(fieldName_, false))
    {
        if (!zoneSubSetPtr_)
        {
            return store
            (
                resultName_,
                fvc::div(lookupObject<FieldType>(fieldName_))
            );
        }
        else
        {
            const fvMeshSubset& subFvMesh = zoneSubSetPtr_->subSetMesh();

            auto tresult
            (
                fvc::div
                (
                    subFvMesh.interpolate
                    (
                        lookupObject<FieldType>(fieldName_),
                        false
                    )
                )
            );

            return storeInDb
            (
                resultName_,
                tresult,
                subFvMesh.subMesh().thisDb()
            );
        }
    }

    return false;
}


template<class Type>
bool Foam::functionObjects::div::writeField()
{
    typedef GeometricField<Type, fvPatchField, volMesh> volField;

    const fvMesh& subFvMesh = zoneSubSetPtr_->subSetMesh().subMesh();
    const auto* fieldPtr = subFvMesh.findObject<volField>(resultName_);

    if (fieldPtr)
    {
        tmp<volField> tmapField = zoneSubSetPtr_->mapToZone<Type>
        (
            subFvMesh.lookupObject<volField>(resultName_)
        );
        tmapField().write();

        return true;
    }

    return false;
}


// ************************************************************************* //
