/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "ChargeParcel.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::ChargeParcel<ParcelType>::propertyList_ =
    Foam::ChargeParcel<ParcelType>::propertyList();


template<class ParcelType>
const std::size_t Foam::ChargeParcel<ParcelType>::sizeofFields
(
    sizeof(ChargeParcel<ParcelType>)
  - offsetof(ChargeParcel<ParcelType>, eq_)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ChargeParcel<ParcelType>::ChargeParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    ParcelType(mesh, is, readFields, newFormat),
    eq_(0.0)
{
    if (readFields)
    {
        if (is.format() == IOstreamOption::ASCII)
        {
            is  >> eq_;
        }
        else if (!is.checkLabelSize<>() || !is.checkScalarSize<>())
        {
            // Non-native label or scalar size

            is.beginRawRead();

            readRawScalar(is, &eq_);

            is.endRawRead();
        }
        else
        {
            is.read(reinterpret_cast<char*>(&eq_), sizeofFields);
        }
    }

    is.check(FUNCTION_NAME);
}


template<class ParcelType>
template<class CloudType>
void Foam::ChargeParcel<ParcelType>::readFields(CloudType& c)
{
    const bool valid = c.size();

    ParcelType::readFields(c);

    IOField<scalar> eq
    (
        c.fieldIOobject("eq", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, eq);

    label i = 0;

    for (ChargeParcel<ParcelType>& p : c)
    {
        p.eq_ = eq[i];

        ++i;
    }
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::ChargeParcel<ParcelType>::readFields
(
    CloudType& c,
    const CompositionType& compModel
)
{
    const bool valid = c.size();

    ParcelType::readFields(c, compModel);

    IOField<scalar> eq
    (
        c.fieldIOobject("eq", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, eq);

    label i = 0;

    for (ChargeParcel<ParcelType>& p : c)
    {
        p.eq_ = eq[i];

        ++i;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::ChargeParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    const label np = c.size();
    const bool valid = np;

    IOField<scalar> eq(c.fieldIOobject("eq", IOobject::NO_READ), np);

    label i = 0;

    for (const ChargeParcel<ParcelType>& p : c)
    {
        eq[i] = p.eq();

        ++i;
    }

    eq.write(valid);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::ChargeParcel<ParcelType>::writeFields
(
    const CloudType& c,
    const CompositionType& compModel
)
{
    ParcelType::writeFields(c, compModel);

    const label np = c.size();
    const bool valid = np;

    IOField<scalar> eq(c.fieldIOobject("eq", IOobject::NO_READ), np);

    label i = 0;

    for (const ChargeParcel<ParcelType>& p : c)
    {
        eq[i] = p.eq();

        ++i;
    }

    eq.write(valid);
}


template<class ParcelType>
void Foam::ChargeParcel<ParcelType>::writeProperties
(
    Ostream& os,
    const wordRes& filters,
    const word& delim,
    const bool namesOnly
) const
{
    ParcelType::writeProperties(os, filters, delim, namesOnly);

    #undef  writeProp
    #define writeProp(Name, Value)                                            \
        ParcelType::writeProperty(os, Name, Value, namesOnly, delim, filters)

    writeProp("eq", eq_);

    #undef writeProp
}


template<class ParcelType>
template<class CloudType>
void Foam::ChargeParcel<ParcelType>::readObjects
(
    CloudType& c,
    const objectRegistry& obr
)
{
    ParcelType::readObjects(c, obr);

    if (!c.size()) return;

    const auto& eq = cloud::lookupIOField<scalar>("eq", obr);

    label i = 0;

    for (ChargeParcel<ParcelType>& p : c)
    {
        p.eq_ = eq[i];

        ++i;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::ChargeParcel<ParcelType>::writeObjects
(
    const CloudType& c,
    objectRegistry& obr
)
{
    ParcelType::writeObjects(c, obr);

    const label np = c.size();

    auto& eq = cloud::createIOField<scalar>("eq", np, obr);

    label i = 0;

    for (const ChargeParcel<ParcelType>& p : c)
    {
        eq[i] = p.eq();

        ++i;
    }
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::ChargeParcel<ParcelType>::readObjects
(
    CloudType& c,
    const CompositionType& compModel,
    const objectRegistry& obr
)
{
    ParcelType::readObjects(c, obr);

    if (!c.size()) return;

    const auto& eq = cloud::lookupIOField<scalar>("eq", obr);

    label i = 0;

    for (ChargeParcel<ParcelType>& p : c)
    {
        p.eq_ = eq[i];

        ++i;
    }
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::ChargeParcel<ParcelType>::writeObjects
(
    const CloudType& c,
    const CompositionType& compModel,
    objectRegistry& obr
)
{
    ParcelType::writeObjects(c, obr);

    const label np = c.size();

    auto& eq = cloud::createIOField<scalar>("eq", np, obr);

    label i = 0;

    for (const ChargeParcel<ParcelType>& p : c)
    {
        eq[i] = p.eq();

        ++i;
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ChargeParcel<ParcelType>& p
)
{
    if (os.format() == IOstreamOption::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.eq();
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.eq_),
            KinematicParcel<ParcelType>::sizeofFields
        );
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
