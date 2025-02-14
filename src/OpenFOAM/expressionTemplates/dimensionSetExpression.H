/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 M. Janssens
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

Namespace
    Foam::dimensionSetExpression

Description
    Expression templates for dimensionSet

SourceFiles
    dimensionSetExpression.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_dimensionSetExpression_H
#define Foam_dimensionSetExpression_H

#include "GenericExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Expression
{

/*---------------------------------------------------------------------------*\
                       Namespace Expression Declarations
\*---------------------------------------------------------------------------*/

// Wrap of non-const reference to dimensionSet
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class dimensionSetRefWrap
:
    public GenericExpression<dimensionSetRefWrap>
{
    dimensionSet& elems_;

public:
    static constexpr bool is_leaf = false;  //true;

    //- Construct from components
    dimensionSetRefWrap(dimensionSet& elems)
    :
        elems_(elems)
    {}

    // Construct from GenericExpression, forcing its evaluation.
    template<typename E>
    dimensionSetRefWrap
    (
        dimensionSet& elems,
        const GenericExpression<E>& expr
    )
    :
        elems_(elems)
    {
        elems_ = expr.evaluate();
    }

    auto evaluate() const
    {
        return elems_;
    }

    template<typename E>
    auto& evaluate
    (
        const GenericExpression<E>& expr
    )
    {
        elems_ = expr.evaluate();
        return elems_;
    }

    //- Assignment
    template<typename E>
    void operator=
    (
        const GenericExpression<E>& expr
    )
    {
        elems_ = expr.evaluate();
    }
};


// Wrap of const reference to dimensionSet
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class dimensionSetConstRefWrap
:
    public GenericExpression<dimensionSetConstRefWrap>
{
    const dimensionSet& elems_;

public:

    // ! Store as copy (since holds reference)
    static constexpr bool is_leaf = false;

    // construct from components
    dimensionSetConstRefWrap(const dimensionSet& elems)
    :
        elems_(elems)
    {}

    auto evaluate() const
    {
        return elems_;
    }
};


} // End namespace Expression

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
