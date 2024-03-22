/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "volFields.H"
#include "volMesh.H"
#include "surfaceFields.H"
#include "surfaceMesh.H"
#include "pointFields.H"
#include "pointMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::PatchFunction1Types::LookupField<Type>::LookupField
(
    const polyPatch& pp,
    const word& redirectType,
    const word& entryName,
    const dictionary& dict,
    const bool faceValues
)
:
    PatchFunction1<Type>(pp, entryName, dict, faceValues),
    lookupBase(dict)
{}


template<class Type>
Foam::PatchFunction1Types::LookupField<Type>::LookupField
(
    const LookupField<Type>& rhs,
    const polyPatch& pp
)
:
    PatchFunction1<Type>(rhs, pp),
    lookupBase(rhs)
{}


template<class Type>
Foam::PatchFunction1Types::LookupField<Type>::LookupField
(
    const LookupField<Type>& rhs
)
:
    LookupField<Type>(rhs, rhs.patch())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::PatchFunction1Types::LookupField<Type>::value(const scalar x) const
{
    const objectRegistry& db = patchFunction1Base::obr();
    const label patchi = patchFunction1Base::patch().index();

    if (this->faceValues())
    {
        typedef GeometricField<Type, fvPatchField, volMesh> VType;
        typedef GeometricField<Type, fvsPatchField, surfaceMesh> SType;

        const auto* fldPtr = db.getObjectPtr<VType>(name_);

        if (fldPtr)
        {
            return fldPtr->boundaryField()[patchi];
        }
        else
        {
            const auto* fldPtr = db.getObjectPtr<SType>(name_);

            if (!fldPtr)
            {
                FatalErrorInFunction
                    << nl
                    << "    failed lookup of " << name_
                    << " (objectRegistry "
                    << db.name()
                    << ")\n    available objects of type " << VType::typeName
                    << ':' << nl
                    << db.names<VType>() << nl
                    << "    available objects of type " << SType::typeName
                    << ':' << nl
                    << db.names<SType>() << nl
                    << exit(FatalError);
            }
            return fldPtr->boundaryField()[patchi];
        }
    }
    else
    {
        // Assume pointField
        typedef GeometricField<Type, pointPatchField, pointMesh> PType;
        const auto& fld = db.lookupObject<PType>(name_);
        return fld.boundaryField()[patchi].patchInternalField();
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::PatchFunction1Types::LookupField<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    return (x2-x1)*value(0.5*(x1+x2));
}


template<class Type>
void Foam::PatchFunction1Types::LookupField<Type>::writeData
(
    Ostream& os
) const
{
    PatchFunction1<Type>::writeData(os);
    os.writeEntry(this->name(), type());
    os.beginBlock(word(this->name() + "Coeffs"));
    writeEntries(os);
    os.endBlock();
}


// ************************************************************************* //
