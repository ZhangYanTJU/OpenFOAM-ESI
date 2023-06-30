/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
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

#include "exprValue.H"
#include "ITstream.H"
#include "Switch.H"
#include <cstring>  // For memcpy, memset

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
static void fillTokens(const Type& val, tokenList& toks)
{
    const direction nCmpt = pTraits<Type>::nComponents;
    const direction nParen = 2*(pTraits<Type>::rank || (nCmpt > 1) ? 1 : 0);

    toks.resize_nocopy(nCmpt + nParen);

    auto iter = toks.begin();

    if (nParen)
    {
        *iter = token::BEGIN_LIST;
        ++iter;
    }

    for (direction cmpt = 0; cmpt < nCmpt; ++cmpt)
    {
        *iter = component(val, cmpt);
        ++iter;
    }

    if (nParen)
    {
        *iter = token::END_LIST;
        ++iter;
    }
}


//- Specialized for bool
template<>
void fillTokens<bool>(const bool& val, tokenList& toks)
{
    toks.resize_nocopy(1);
    toks.front() = token::boolean(val);
}


//- Specialized for label
template<>
void fillTokens<label>(const label& val, tokenList& toks)
{
    toks.resize_nocopy(1);
    toks.front() = val;
}


//- Specialized for scalar
template<>
void fillTokens<scalar>(const scalar& val, tokenList& toks)
{
    toks.resize_nocopy(1);
    toks.front() = val;
}

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::expressions::valueTypeCode
Foam::expressions::exprValue::peekType(const ITstream& is)
{
    expressions::valueTypeCode whichCode(expressions::valueTypeCode::INVALID);

    const token& tok0 = is.peek();

    if (tok0.isPunctuation(token::BEGIN_LIST))
    {
        // Expecting "( content )" - eg, (x y z), (xx xy ...)

        // First component starts after the opening '('.
        // Can use the current index if the '(' actually came from
        // the putBack.

        const label firstCmpti = (is.tokenIndex() + (is.hasPutback() ? 0 : 1));

        // Search for closing ')', require all components to be numbers
        for (label endCmpti = firstCmpti; endCmpti < is.size(); ++endCmpti)
        {
            const token& tok = is[endCmpti];

            ///InfoErr
            ///    << "check token: " << (endCmpti-firstCmpti) << " : "
            ///    << is[endCmpti].name() << nl;

            if (tok.isPunctuation(token::END_LIST))
            {
                // Select based on the number of components
                // cf. pTraits<Type>::nComponents

                switch (endCmpti - firstCmpti)
                {
                    case 1:  // pTraits<sphericalTensor>::nComponents
                        whichCode = exprTypeTraits<sphericalTensor>::value;
                        break;

                    // FUTURE?
                    // case 2:  // pTraits<complex>::nComponents
                    //     whichCode = exprTypeTraits<complex>::value;
                    //     break;

                    case 3:  // pTraits<vector>::nComponents
                        whichCode = exprTypeTraits<vector>::value;
                        break;

                    case 6:  // pTraits<symmTensor>::nComponents
                        whichCode = exprTypeTraits<symmTensor>::value;
                        break;

                    case 9:  // pTraits<tensor>::nComponents
                        whichCode = exprTypeTraits<tensor>::value;
                        break;
                }

                // Closing ')' terminates peeking
                break;
            }
            else if (!tok.isNumber())
            {
                // All components should be numeric
                break;
            }
        }
    }
    else if (tok0.good())
    {
        /// InfoErr<< "check token: " << tok0.info() << nl;

        if (tok0.isScalar())
        {
            whichCode = exprTypeTraits<scalar>::value;
        }
        else if (tok0.isLabel())
        {
            whichCode = exprTypeTraits<label>::value;
        }
        else if (Switch(tok0).good())
        {
            whichCode = exprTypeTraits<bool>::value;
        }
    }

    return whichCode;
}


bool Foam::expressions::exprValue::read
(
    const std::string& str,
    exprValue& val
)
{
    ITstream is(str);

    // No trailing non-whitespace!
    return (val.read(is) && !is.nRemainingTokens());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::expressions::exprValue::clear()
{
    std::memset(static_cast<void*>(this), '\0', sizeof(*this));
    // Redundant:  typeCode_ = expressions::valueTypeCode::NONE;
}


void Foam::expressions::exprValue::deepCopy(const exprValue& rhs)
{
    if (this != &rhs)
    {
        // Self-assignment is a no-op
        std::memcpy(static_cast<void*>(this), &rhs, sizeof(*this));
    }
}


Foam::tokenList Foam::expressions::exprValue::tokens() const
{
    tokenList toks;

    switch (typeCode_)
    {
        #undef  doLocalCode
        #define doLocalCode(Type, UnusedParam)                        \
                                                                      \
        case expressions::valueTypeCode::type_##Type :                \
        {                                                             \
            const Type* dataPtr = data_.get<Type>();                  \
            if (dataPtr)                                              \
            {                                                         \
                fillTokens<Type>(*dataPtr, toks);                     \
            }                                                         \
            break;                                                    \
        }

        FOR_ALL_EXPR_VALUE_TYPES(doLocalCode);
        #undef doLocalCode

        // exprValue may only be a subset of valueTypeCode types
        default: break;
    }

    return toks;
}


void Foam::expressions::exprValue::print(Ostream& os) const
{
    switch (typeCode_)
    {
        #undef  doLocalCode
        #define doLocalCode(Type, UnusedParam)                        \
                                                                      \
        case expressions::valueTypeCode::type_##Type :                \
        {                                                             \
            const Type* dataPtr = data_.get<Type>();                  \
            if (dataPtr)                                              \
            {                                                         \
                os << *dataPtr;                                       \
            }                                                         \
            break;                                                    \
        }

        FOR_ALL_EXPR_VALUE_TYPES(doLocalCode);
        #undef doLocalCode

        // exprValue may only be a subset of valueTypeCode types
        default: break;
    }
}


bool Foam::expressions::exprValue::read(ITstream& is)
{
    clear();  // type: none, value: zero

    const valueTypeCode whichCode(exprValue::peekType(is));

    switch (whichCode)
    {
        #undef  doLocalCode
        #define doLocalCode(Type, UnusedParam)                        \
                                                                      \
        case expressions::valueTypeCode::type_##Type :                \
        {                                                             \
            data_.set<Type>(pTraits<Type>(is));                       \
            typeCode_ = whichCode;                                    \
            return true;                                              \
        }

        FOR_ALL_EXPR_VALUE_TYPES(doLocalCode);
        #undef doLocalCode

        // exprValue may only be a subset of valueTypeCode types
        default: break;
    }

    return false;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const expressions::exprValue& val
)
{
    val.print(os);
    return os;
}


template<>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<expressions::exprValue>& iproxy
)
{
    const auto& val = *iproxy;

    if (val.good())
    {
        os  << val.valueTypeName() << ": ";
        val.print(os);
    }
    else
    {
        // typeCode_ is *never* set to INVALID,
        // so NONE is the only remaining non-good type
        os  << "none";
    }

    return os;
}


// ************************************************************************* //
