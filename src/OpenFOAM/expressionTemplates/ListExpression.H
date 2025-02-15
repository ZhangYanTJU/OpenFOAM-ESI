/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024,2025 M. Janssens
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
    Foam::ListExpression

Description
    Expression templates for List

SourceFiles
    ListExpression.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_ListExpression_H
#define Foam_ListExpression_H

#include <cassert>
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Expression
{

/*---------------------------------------------------------------------------*\
                       Namespace Expression Declarations
\*---------------------------------------------------------------------------*/

template<typename E, typename ValType>
class ListExpression
{
public:
    static constexpr bool is_leaf = false;

    ValType operator[](const label i) const
    {
        // Delegation to the actual expression type. This avoids dynamic
        // polymorphism (a.k.a. virtual functions in C++)
        return static_cast<E const&>(*this)[i];
    }
    label size() const { return static_cast<E const&>(*this).size(); }

    //- Helper: assign to passed in list
    template<class Container>
    Container& evaluate(Container& lst) const
    {
        //Pout<< "in ListExpression evaluate with Container:" << lst.size()
        //    << endl;

        lst.resize_nocopy(this->size());
        for (label i = 0; i < lst.size(); ++i)
        {
            lst[i] = operator[](i);
        }
        return lst;
    }
};


// Wrap of List
// ~~~~~~~~~~~~

template<class T>
class ListWrap
:
    public ListExpression<ListWrap<T>, T>
{
    List<T> elems_;

public:
    static constexpr bool is_leaf = true;

    typedef T value_type;

    //- Copy construct
    ListWrap(const ListWrap<T>&) = default;

    //- Move construct
    ListWrap(ListWrap<T>&&) = default;

    // Construct from initializer list 
    ListWrap(std::initializer_list<T> init)
    :
        elems_(init)
    {}

    // Construct from ListExpression, forcing its evaluation.
    template<typename E>
    ListWrap(const ListExpression<E, typename E::value_type>& expr)
    :
        elems_(expr.size())
    {
        expr.evaluate(elems_);
    }

    // Access underlying data
    const auto& data() const
    { 
        return elems_; 
    }

    // Access underlying data
    auto& data()
    { 
        return elems_; 
    }

    value_type operator[](const label i) const
    {
        return elems_[i];
    }

    value_type& operator[](const label i)
    {
        return elems_[i];
    }

    label size() const
    {
        return elems_.size();
    }
};


// Wrap of non-const reference to List
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template<class T>
class ListRefWrap
:
    public ListExpression<ListRefWrap<T>, T>
{
    List<T>& elems_;

public:
    static constexpr bool is_leaf = false;  //true;

    typedef T value_type;

    //- Copy construct
    ListRefWrap(ListRefWrap<T>& w)
    :
        elems_(w.elems_)
    {
        //Pout<< "    constructed copy listRefWrap for size:" << elems_.size()
        //    << endl;
    }

    //- Move construct
    ListRefWrap(ListRefWrap<T>&& elems)
    :
        elems_(elems)
    {}

    //- Construct from List and expected size. This is to enforce that
    //- our size() member function returns the correct size. Use the
    //- construct-from-expression otherwise.
    ListRefWrap(const label size, List<T>& elems)
    :
        elems_(elems)
    {
        elems_.resize_nocopy(size);

        //Pout<< "    ListRefWrap : constructed from List:"
        //    << " listRefWrap for size:" << elems_.size() << endl;
    }

    // Construct from ListExpression, forcing its evaluation.
    template<typename E>
    ListRefWrap
    (
        List<T>& elems,
        const ListExpression<E, typename E::value_type>& expr
    )
    :
        elems_(elems)
    {
        //Pout<< "    ListRefWrap : constructed from expression size:"
        //    << expr.size() << endl;

        expr.evaluate(elems_);

    }

    //- Assignment
    template<typename E>
    void operator=(const ListExpression<E, typename E::value_type>& expr)
    {
        expr.evaluate(elems_);
    }

    //- Evaluate and return as GeoField. Rename to evaluate to make it clear
    //- it takes time? Or leave as indexing for convenience?
    template<typename E>
    List<T>& evaluate(const ListExpression<E, typename E::value_type>& expr)
    {
        return expr.evaluate(elems_);
    }

    // Access underlying data
    const auto& data() const
    { 
        return elems_; 
    }

    // Access underlying data
    auto& data()
    { 
        return elems_; 
    }

    value_type operator[](const label i) const
    {
        return elems_[i];
    }

    value_type& operator[](const label i)
    {
        return elems_[i];
    }

    label size() const
    {
        return elems_.size();
    }
};


// Wrap of const reference to List
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template<class T>
class ListConstRefWrap
:
    public ListExpression<ListConstRefWrap<T>, T>
{
    const List<T>& elems_;

public:

    // ! Store as copy (since holds reference)
    static constexpr bool is_leaf = false;

    typedef T value_type;

    // construct from components
    ListConstRefWrap(const List<T>& elems)
    :
        elems_(elems)
    {
        //Pout<< "    ListConstRefWrap :"
        //    << " constructed with elems_.size:" << elems_.size() << endl;
    }

    // returns the underlying data
    const List<T>& data() const
    { 
        return elems_; 
    }

    value_type operator[](const label i) const
    {
        return elems_[i];
    }
    label size() const
    {
        return elems_.size();
    }
};


// Expresssions on List
// ~~~~~~~~~~~~~~~~~~~~

#define EXPRESSION_FUNCTION1(Func, BaseFunc, OpFunc)                           \
template<typename E1>                                                          \
class OpFunc                                                                   \
:                                                                              \
    public ListExpression<OpFunc<E1>, typename E1::value_type>                 \
{                                                                              \
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type u_;      \
                                                                               \
public:                                                                        \
    static constexpr bool is_leaf = false;                                     \
                                                                               \
    typedef typename E1::value_type value_type;                                \
                                                                               \
    OpFunc(const E1& u)                                                        \
    :                                                                          \
        u_(u)                                                                  \
    {}                                                                         \
                                                                               \
    decltype(auto) operator[](const label i) const                             \
    {                                                                          \
        return BaseFunc(u_[i]);                                                \
    }                                                                          \
                                                                               \
    label size() const { return u_.size(); }                                   \
};                                                                             \
template<typename E1>                                                          \
OpFunc<E1> Func                                                                \
(                                                                              \
    const ListExpression<E1, typename E1::value_type>& u                       \
)                                                                              \
{                                                                              \
    return OpFunc<E1>(static_cast<const E1&>(u));                              \
}


// Do '-' separately until we work out macro expansion...
template<typename E1>
class List_negate
:
    public ListExpression<List_negate<E1>, typename E1::value_type>
{
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type u_;

public:
    static constexpr bool is_leaf = false;

    typedef typename E1::value_type value_type;

    List_negate(const E1& u)
    :
        u_(u)
    {}

    decltype(auto) operator[](const label i) const
    {
        return -u_[i];
    }

    label size() const { return u_.size(); }
};
template<typename E1>
List_negate<E1> operator-
(
    const ListExpression<E1, typename E1::value_type>& u
)
{
    return List_negate<E1>(static_cast<const E1&>(u));
}


#define EXPRESSION_FUNCTION2(Func, BaseFunc, OpFunc)                           \
template<typename E1, typename E2>                                             \
class OpFunc                                                                   \
:                                                                              \
    public ListExpression                                                      \
    <                                                                          \
        OpFunc<E1, E2>,                                                        \
        typename E1::value_type                                                \
    >                                                                          \
{                                                                              \
    /* cref if leaf, copy otherwise */                                         \
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type u_;      \
    typename std::conditional<E2::is_leaf, const E2&, const E2>::type v_;      \
                                                                               \
public:                                                                        \
    static constexpr bool is_leaf = false;                                     \
                                                                               \
    /* Type to return for internal field */                                    \
    typedef typename E1::value_type value_type;                                \
                                                                               \
    OpFunc(const E1& u, const E2& v)                                           \
    :                                                                          \
     u_(u), v_(v)                                                              \
    {                                                                          \
        assert(u.size() == v.size());                                          \
    }                                                                          \
    decltype(auto) operator[](const label i) const                             \
    {                                                                          \
        return BaseFunc(u_[i], v_[i]);                                         \
    }                                                                          \
    label size() const {return u_.size(); }                                    \
};                                                                             \
template<typename E1, typename E2>                                             \
OpFunc<E1, E2> Func                                                            \
(                                                                              \
    ListExpression                                                             \
    <                                                                          \
        E1,                                                                    \
        typename E1::value_type                                                \
    > const& u,                                                                \
    ListExpression                                                             \
    <                                                                          \
        E2,                                                                    \
        typename E2::value_type                                                \
    > const& v                                                                 \
)                                                                              \
{                                                                              \
    return OpFunc<E1, E2>                                                      \
    (                                                                          \
        static_cast<const E1&>(u),                                             \
        static_cast<const E2&>(v)                                              \
    );                                                                         \
}

#define EXPRESSION_OPERATOR(Op, OpFunc)                                        \
template<typename E1, typename E2>                                             \
class OpFunc                                                                   \
:                                                                              \
    public ListExpression                                                      \
    <                                                                          \
        OpFunc<E1, E2>,                                                        \
        typename E1::value_type                                                \
    >                                                                          \
{                                                                              \
    /* cref if leaf, copy otherwise */                                         \
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type u_;      \
    typename std::conditional<E2::is_leaf, const E2&, const E2>::type v_;      \
                                                                               \
public:                                                                        \
    static constexpr bool is_leaf = false;                                     \
                                                                               \
    /* Type to return for internal field */                                    \
    typedef typename E1::value_type value_type;                                \
                                                                               \
    OpFunc(const E1& u, const E2& v)                                           \
    :                                                                          \
     u_(u), v_(v)                                                              \
    {                                                                          \
        assert(u.size() == v.size());                                          \
    }                                                                          \
    decltype(auto) operator[](const label i) const                             \
    {                                                                          \
        return u_[i] Op v_[i];                                                 \
    }                                                                          \
    label size() const { return u_.size(); }                                   \
};                                                                             \
template<typename E1, typename E2>                                             \
OpFunc<E1, E2>                                                                 \
operator Op                                                                    \
(                                                                              \
    ListExpression                                                             \
    <                                                                          \
        E1,                                                                    \
        typename E1::value_type                                                \
    > const& u,                                                                \
    ListExpression                                                             \
    <                                                                          \
        E2,                                                                    \
        typename E2::value_type                                                \
    > const& v                                                                 \
)                                                                              \
{                                                                              \
    return OpFunc<E1, E2>                                                      \
    (                                                                          \
        static_cast<const E1&>(u),                                             \
        static_cast<const E2&>(v)                                              \
    );                                                                         \
}


//- List functions
// macro arguments:
//  - 'sin' : name of new function (in ListExpression namespace)
//  - '::sin' : per-element function to call
//  - 'List_sin' : generated expression helper class
EXPRESSION_FUNCTION1(sin, ::sin, List_sin)
EXPRESSION_FUNCTION1(cos, ::cos, List_cos)
EXPRESSION_FUNCTION1(tan, ::tan, List_tan)
EXPRESSION_FUNCTION1(sinh, ::sinh, List_sinh)
EXPRESSION_FUNCTION1(cosh, ::cosh, List_cosh)
EXPRESSION_FUNCTION1(tanh, ::tanh, List_tanh)
EXPRESSION_FUNCTION1(sqr, Foam::sqr, List_sqr)
EXPRESSION_FUNCTION1(sqrt, Foam::sqrt, List_sqrt)
EXPRESSION_FUNCTION1(magSqr, Foam::magSqr, List_magSqr)
EXPRESSION_FUNCTION1(symm, Foam::symm, List_symm)
//EXPRESSION_FUNCTION1(-, -, List_negate)
EXPRESSION_FUNCTION1(pow2, Foam::pow2, List_pow2)
EXPRESSION_FUNCTION1(pow3, Foam::pow3, List_pow3)
EXPRESSION_FUNCTION1(pow4, Foam::pow4, List_pow4)

#undef EXPRESSION_FUNCTION1

EXPRESSION_FUNCTION2(min, Foam::min, List_min)
EXPRESSION_FUNCTION2(max, Foam::max, List_max)

#undef EXPRESSION_FUNCTION2

//- List operators
EXPRESSION_OPERATOR(+, List_add)
EXPRESSION_OPERATOR(-, List_subtract)
EXPRESSION_OPERATOR(*, List_multiply)
EXPRESSION_OPERATOR(/, List_divide)

#undef EXPRESSION_OPERATOR

//- Not possible: a *= b has to return a new type so it can never be in-place
//EXPRESSION_LIST_OPERATOR(+=, List_inplace_add)
//EXPRESSION_LIST_OPERATOR(-=, List_inplace_subtract)
//EXPRESSION_LIST_OPERATOR(*=, List_inplace_multiply)
//EXPRESSION_LIST_OPERATOR(/=, List_inplace_divide)
//#undef EXPRESSION_LIST_OPERATOR


// Expressions on constants
// ~~~~~~~~~~~~~~~~~~~~~~~~

template<class T>
class constantListWrap
:
    public ListExpression<constantListWrap<T>, T>
{
private:

    const label size_;

    const T val_;

public:

    static constexpr bool is_leaf = false;

    typedef T value_type;

    // construct from components
    constantListWrap(const label size, const T val)
    :
        size_(size),
        val_(val)
    {
        //Pout<< "    constantListWrap : constructed from size:"
        //    << size_ << " value:" << val_ << endl;
    }

    // construct from components
    constantListWrap(const label size, const dimensioned<T>& val)
    :
        size_(size),
        val_(val.value())
    {
        //Pout<< "    constantListWrap : constructed from size:"
        //    << size_ << " dimensionedvalue:" << val_ << endl;
    }

    T operator[](const label i) const
    {
        return val_;
    }
    label size() const
    {
        //Pout<< "    constantListWrap : returning size:"
        //    << size_ << " for dimensionedvalue:" << val_ << endl;
        return size_;
    }
};


} // End namespace Expression

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//    #include "ListExpression.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
