/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2025 OpenCFD Ltd.
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

#include "DimensionedField.H"
#include "IOstreams.H"
#include "localIOdictionary.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, class GeoMesh>
void Foam::DimensionedField<Type, GeoMesh>::readField
(
    const dictionary& fieldDict,
    const word& fieldDictEntry
)
{
    dimensions_.readEntry("dimensions", fieldDict);

    // Note: oriented state may have already been set on construction
    // - if so - do not reset by re-reading
    // - required for backwards compatibility in case of restarting from
    //   an old run where the oriented state may not have been set
    if (oriented_.oriented() != orientedType::ORIENTED)
    {
        oriented_.read(fieldDict);  // The "oriented" entry (if present)
    }


    // The primitive field
    auto& fld = static_cast<Field<Type>&>(*this);

    fld.resize_nocopy(GeoMesh::size(mesh_));
    fld.assign(fieldDictEntry, fieldDict, fld.size());  // <- MUST_READ
}


template<class Type, class GeoMesh>
void Foam::DimensionedField<Type, GeoMesh>::readField
(
    const word& fieldDictEntry
)
{
    dictionary dict
    (
        localIOdictionary::readContents
        (
            IOobject
            (
                this->name(),
                this->instance(),
                this->local(),
                this->db(),
                IOobjectOption::MUST_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::NO_REGISTER
            ),
            typeName
        )
    );

    this->close();

    readField(dict, fieldDictEntry);
}


template<class Type, class GeoMesh>
bool Foam::DimensionedField<Type, GeoMesh>::readIfPresent
(
    const word& fieldDictEntry
)
{
    if
    (
        this->isReadRequired()
     || (this->isReadOptional() && this->headerOk())
    )
    {
        readField(fieldDictEntry);
        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Foam::DimensionedField<Type, GeoMesh>::DimensionedField
(
    const IOobject& io,
    const Mesh& mesh,
    const word& fieldDictEntry
)
:
    regIOobject(io),
    mesh_(mesh)
{
    readField(fieldDictEntry);
}


template<class Type, class GeoMesh>
Foam::DimensionedField<Type, GeoMesh>::DimensionedField
(
    const IOobject& io,
    const Mesh& mesh,
    const dictionary& fieldDict,
    const word& fieldDictEntry
)
:
    regIOobject(io),
    mesh_(mesh)
{
    readField(fieldDict, fieldDictEntry);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
bool Foam::DimensionedField<Type, GeoMesh>::writeData
(
    Ostream& os,
    const word& fieldDictEntry
) const
{
    os.writeEntry("dimensions", dimensions());
    os << nl;

    if (oriented_.writeEntry(os))
    {
        os << nl;
    }

    Field<Type>::writeEntry(fieldDictEntry, os);

    os.check(FUNCTION_NAME);
    return os.good();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const DimensionedField<Type, GeoMesh>& fld
)
{
    fld.writeData(os);

    return os;
}


template<class Type, class GeoMesh>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const tmp<DimensionedField<Type, GeoMesh>>& tfld
)
{
    tfld().writeData(os);
    tfld.clear();

    return os;
}


// ************************************************************************* //
