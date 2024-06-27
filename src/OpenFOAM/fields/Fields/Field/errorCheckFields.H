/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

Description
    Field size checks for field operations

\*---------------------------------------------------------------------------*/

#ifndef Foam_errorCheckFields_H
#define Foam_errorCheckFields_H

#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type1, class Type2>
void checkFields
(
    const UList<Type1>& f1,
    const UList<Type2>& f2,
    const char* op
)
{
    #ifdef FULLDEBUG
    if (f1.size() != f2.size())
    {
        FatalErrorInFunction
            << " Field<"<<pTraits<Type1>::typeName<<"> f1("<<f1.size()<<')'
            << " and Field<"<<pTraits<Type2>::typeName<<"> f2("<<f2.size()<<')'
            << nl
            << " for operation " << op << endl
            << abort(FatalError);
    }
    #endif
}


template<class Type1, class Type2, class Type3>
void checkFields
(
    const UList<Type1>& f1,
    const UList<Type2>& f2,
    const UList<Type3>& f3,
    const char* op
)
{
    #ifdef FULLDEBUG
    if
    (
        f1.size() != f2.size()
     || f1.size() != f3.size()
    )
    {
        FatalErrorInFunction
            << " Field<"<<pTraits<Type1>::typeName<<"> f1("<<f1.size()<<')'
            << ", Field<"<<pTraits<Type2>::typeName<<"> f2("<<f2.size()<<')'
            << " and Field<"<<pTraits<Type3>::typeName<<"> f3("<<f3.size()<<')'
            << nl
            << " for operation " << op << endl
            << abort(FatalError);
    }
    #endif
}


template<class Type1, class Type2, class Type3, class Type4>
void checkFields
(
    const UList<Type1>& f1,
    const UList<Type2>& f2,
    const UList<Type3>& f3,
    const UList<Type4>& f4,
    const char* op
)
{
    #ifdef FULLDEBUG
    if
    (
        f1.size() != f2.size()
     || f1.size() != f3.size()
     || f1.size() != f4.size()
    )
    {
        FatalErrorInFunction
            << " Field<"<<pTraits<Type1>::typeName<<"> f1("<<f1.size()<<')'
            << ", Field<"<<pTraits<Type2>::typeName<<"> f2("<<f2.size()<<')'
            << ", Field<"<<pTraits<Type3>::typeName<<"> f3("<<f3.size()<<')'
            << " and Field<"<<pTraits<Type4>::typeName<<"> f4("<<f4.size()<<')'
            << nl
            << " for operation " << op << endl
            << abort(FatalError);
    }
    #endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
