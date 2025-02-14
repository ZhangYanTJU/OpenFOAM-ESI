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
    Foam::GenericExpression

Description
    Expression templates.

SourceFiles
    GenericExpression.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_GenericExpression_H
#define Foam_GenericExpression_H

#include <cassert>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Expression
{

/*---------------------------------------------------------------------------*\
                     Class GenericExpression Declarations
\*---------------------------------------------------------------------------*/

template<typename E>
class GenericExpression
{
public:
    static constexpr bool is_leaf = false;

    auto evaluate() const
    {
        // Delegation to the actual expression type. This avoids dynamic
        // polymorphism (a.k.a. virtual functions in C++)
        return static_cast<E const&>(*this).evaluate();
    }
};


// Macros
// ~~~~~~

#define EXPRESSION_FUNCTION1(BaseClass, Func, BaseFunc, OpFunc)                \
template<typename E1>                                                          \
class OpFunc                                                                   \
:                                                                              \
    public BaseClass<OpFunc<E1>>                                               \
{                                                                              \
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type u_;      \
                                                                               \
public:                                                                        \
    static constexpr bool is_leaf = false;                                     \
                                                                               \
    OpFunc(const E1& u)                                                        \
    :                                                                          \
        u_(u)                                                                  \
    {}                                                                         \
                                                                               \
    auto evaluate() const                                                      \
    {                                                                          \
        return BaseFunc(u_.evaluate());                                        \
    }                                                                          \
};                                                                             \
template<typename E1>                                                          \
OpFunc<E1> Func                                                                \
(                                                                              \
    const BaseClass<E1>& u                                                     \
)                                                                              \
{                                                                              \
    return OpFunc<E1>(static_cast<const E1&>(u));                              \
}


#define EXPRESSION_FUNCTION2(BaseClass, Func, BaseFunc, OpFunc)                \
template<typename E1, typename E2>                                             \
class OpFunc                                                                   \
:                                                                              \
    public BaseClass<OpFunc<E1, E2>>                                           \
{                                                                              \
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type u_;      \
    typename std::conditional<E2::is_leaf, const E2&, const E2>::type v_;      \
                                                                               \
public:                                                                        \
    static constexpr bool is_leaf = false;                                     \
                                                                               \
    OpFunc(const E1& u, const E2& v)                                           \
    :                                                                          \
     u_(u), v_(v)                                                              \
    {}                                                                         \
    auto evaluate() const                                                      \
    {                                                                          \
        return BaseFunc(u_.evaluate(), v_.evaluate());                         \
    }                                                                          \
};                                                                             \
template<typename E1, typename E2>                                             \
OpFunc<E1, E2> Func                                                            \
(                                                                              \
    const BaseClass<E1>& u,                                                    \
    const BaseClass<E2>& v                                                     \
)                                                                              \
{                                                                              \
    return OpFunc<E1, E2>                                                      \
    (                                                                          \
        static_cast<const E1&>(u),                                             \
        static_cast<const E2&>(v)                                              \
    );                                                                         \
}


#define EXPRESSION_OPERATOR(BaseClass, Op, BaseOp, OpFunc)                     \
template<typename E1, typename E2>                                             \
class OpFunc                                                                   \
:                                                                              \
    public BaseClass                                                           \
    <                                                                          \
        OpFunc<E1, E2>                                                         \
    >                                                                          \
{                                                                              \
    /* cref if leaf, copy otherwise */                                         \
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type u_;      \
    typename std::conditional<E2::is_leaf, const E2&, const E2>::type v_;      \
                                                                               \
public:                                                                        \
    static constexpr bool is_leaf = false;                                     \
                                                                               \
    OpFunc(const E1& u, const E2& v)                                           \
    :                                                                          \
     u_(u), v_(v)                                                              \
    {}                                                                         \
    auto evaluate() const                                                      \
    {                                                                          \
        return u_.evaluate() BaseOp v_.evaluate();                             \
    }                                                                          \
};                                                                             \
template<typename E1, typename E2>                                             \
OpFunc<E1, E2>                                                                 \
operator Op                                                                    \
(                                                                              \
    BaseClass<E1> const& u,                                                    \
    BaseClass<E2> const& v                                                     \
)                                                                              \
{                                                                              \
    return OpFunc<E1, E2>                                                      \
    (                                                                          \
        static_cast<const E1&>(u),                                             \
        static_cast<const E2&>(v)                                              \
    );                                                                         \
}


// Do '-' separately until we work out macro expansion...
template<typename E1>
class g_negate
:
    public GenericExpression
    <
        g_negate<E1>
    >
{
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type u_;

public:
    static constexpr bool is_leaf = false;

    g_negate(const E1& u)
    :
        u_(u)
    {}

    auto evaluate() const
    {
        return -u_.evaluate();
    }
};
template<typename E1>
g_negate<E1> operator-
(
    const GenericExpression<E1>& u
)
{
    return g_negate<E1>(static_cast<const E1&>(u));
}


EXPRESSION_FUNCTION1(GenericExpression, sqr, Foam::sqr, g_sqr)
EXPRESSION_FUNCTION1(GenericExpression, sqrt, Foam::sqrt, g_sqrt)
EXPRESSION_FUNCTION1(GenericExpression, magSqr, Foam::magSqr, g_magSqr)
EXPRESSION_FUNCTION1(GenericExpression, symm, Foam::symm, g_symm)
EXPRESSION_FUNCTION1(GenericExpression, pow2, Foam::pow2, g_pow2)
EXPRESSION_FUNCTION1(GenericExpression, pow3, Foam::pow3, g_pow3)
EXPRESSION_FUNCTION1(GenericExpression, pow4, Foam::pow4, g_pow4)
EXPRESSION_FUNCTION1(GenericExpression, operator~, Foam::operator~, g_compl)

//TBD. Parse problem
//EXPRESSION_FUNCTION1(GenericExpression, operator!, Foam::operator!, g_not)

#undef EXPRESSION_FUNCTION1

EXPRESSION_FUNCTION2(GenericExpression, operator+, Foam::operator+, g_add);
EXPRESSION_FUNCTION2(GenericExpression, operator-, Foam::operator-, g_subtract);
EXPRESSION_FUNCTION2(GenericExpression, operator*, Foam::operator*, g_multiply);
EXPRESSION_FUNCTION2(GenericExpression, operator/, Foam::operator/, g_divide);

#undef EXPRESSION_FUNCTION2

EXPRESSION_OPERATOR(GenericExpression, ||, ||, g_or)
EXPRESSION_OPERATOR(GenericExpression, &&, &&, g_and)
EXPRESSION_OPERATOR(GenericExpression, &, &, g_bitand)
EXPRESSION_OPERATOR(GenericExpression, |, |, g_bitor)
EXPRESSION_OPERATOR(GenericExpression, ^, ^, g_xor)

#undef EXPRESSION_OPERATOR

} // End namespace Expression

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
