/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024,2025 M. Janssens
    Copyright (C) 2025 OpenCFD Ltd.
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
    Foam::Expression

Description
    Expression templates for GeometricFields

SourceFiles
    GeometricFieldExpression.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_GeometricFieldExpression_H
#define Foam_GeometricFieldExpression_H

#include "ListExpression.H"
#include "pointMesh.H"
//#include "fixedValuePointPatchField.H"  // no assignable yet. tbd.
//#include <boost/core/demangle.hpp>
//#include <typeinfo>
//#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Expression
{

/*---------------------------------------------------------------------------*\
                   Namespace Expression Declarations
\*---------------------------------------------------------------------------*/


// Expression of GeometricField
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template
<
    typename E,
    typename IntExpr,
    typename UncoupledPatchExpr,
    typename CoupledPatchExpr,
    typename Type
>
class GeometricFieldExpression
{
protected:

    const dimensionSet dimensions_;
    const orientedType oriented_;

public:
    static constexpr bool is_leaf = false;

    GeometricFieldExpression
    (
        const dimensionSet& dimensions,
        const orientedType oriented
    )
    :
        dimensions_(dimensions),
        oriented_(oriented)
    {
        //Pout<< "GeometricFieldExpression :" << nl
        //    << "    E          : "
        //    << boost::core::demangle(typeid(E).name())
        //    << nl
        //    << "    IntExpr   : "
        //    << boost::core::demangle(typeid(IntExpr).name()) << nl
        //    << "    dimensions : " << dimensions_ << nl
        //    << "    oriented   : " << oriented_ << nl
        //    << endl;
    }

    // TBD. Do we still need this? Is faster than always going through
    // internalField() though.
    Type operator[](const label i) const
    {
        // Delegation to the actual expression type. This avoids dynamic
        // polymorphism (a.k.a. virtual functions in C++)
        return static_cast<E const&>(*this)[i];
    }
    auto size() const noexcept { return static_cast<const E&>(*this).size(); }

    const dimensionSet& dimensions() const noexcept
    {
        return dimensions_;
    }

    const orientedType& oriented() const noexcept
    {
        return oriented_;
    }

    IntExpr internalField() const
    {
        return static_cast<const E&>(*this).internalField();
    }

    IntExpr internalField()
    {
        return static_cast<E&>(*this).internalField();
    }

    UncoupledPatchExpr patchField(const label i)
    {
        return static_cast<E&>(*this).patchField(i);
    }

    UncoupledPatchExpr patchField(const label i) const
    {
        return static_cast<const E&>(*this).patchField(i);
    }

    CoupledPatchExpr coupledPatchField(const label i)
    {
        return static_cast<E&>(*this).coupledPatchField(i);
    }

    CoupledPatchExpr coupledPatchField(const label i) const
    {
        return static_cast<const E&>(*this).coupledPatchField(i);
    }

    template<class Op>
    auto access(const Op& cop, const label i) const
    {
        return static_cast<const E&>(*this).access(cop, i);
    }

    //- Helper to evaluate a GeometricField
    template<class GeoField>
    GeoField& evaluate(GeoField& fld, const bool force = false) const
    {
        fld.dimensions() = this->dimensions();
        fld.oriented() = this->oriented();

        assert(fld.size() == this->size());
        using ActualType = typename IntExpr::value_type;
        ListRefWrap<ActualType> wfld(fld.internalFieldRef(), internalField());

        // Do boundary field
        auto& bfld = fld.boundaryFieldRef();
        const label n = bfld.size();
        for (label i = 0; i < n; ++i)
        {
            auto& pfld = bfld[i];
            if (force || pfld.assignable())
            {
                if (pfld.coupled())
                {
                    // Get patch expression
                    const auto patchExpr = this->coupledPatchField(i);

                    // Note: check can be removed once assignable also on
                    //       pointPatchFields. TBD.
                    if constexpr
                    (
                        std::is_same_v<pointMesh, typename GeoField::Mesh>
                    )
                    {
                        const auto* vp = isA<Field<Type>>(pfld);
                        if (vp)
                        {
                            patchExpr.evaluate(const_cast<Field<Type>&>(*vp));
                        }
                    }
                    else
                    {
                        patchExpr.evaluate(pfld);
                    }
                }
                else
                {
                    // Get patch expression
                    const auto patchExpr = this->patchField(i);

                    // Note: check can be removed once assignable also on
                    //       pointPatchFields. TBD.
                    if constexpr
                    (
                        std::is_same_v<pointMesh, typename GeoField::Mesh>
                    )
                    {
                        const auto* vp = isA<Field<Type>>(pfld);
                        if (vp)
                        {
                            patchExpr.evaluate(const_cast<Field<Type>&>(*vp));
                        }
                    }
                    else
                    {
                        patchExpr.evaluate(pfld);
                    }
                }
            }
        }
        fld.correctLocalBoundaryConditions();
        return fld;
    }
};


template<class GeoField>
class GeometricFieldRefWrap
:
    public GeometricFieldExpression
    <
        GeometricFieldRefWrap<GeoField>,
        ListRefWrap<typename GeoField::value_type>,
        ListRefWrap<typename GeoField::value_type>,
        ListTmpWrap<Field<typename GeoField::value_type>>,
        typename GeoField::value_type
    >
{
public:

    static constexpr bool is_leaf = false;  //true;

    //- The GeometricField type
    typedef GeoField this_type;

    //- Type to return for internal field
    typedef typename GeoField::value_type value_type;

    //- Type to return for patchField
    typedef ListRefWrap<value_type> IntExpr;
    typedef ListRefWrap<value_type> UncoupledPatchExpr;
    typedef ListTmpWrap<Field<value_type>> CoupledPatchExpr;


private:

    this_type& elems_;


public:

    //- Copy construct
    GeometricFieldRefWrap(this_type& elems)
    :
        GeometricFieldExpression
        <
            GeometricFieldRefWrap<GeoField>,
            IntExpr,
            UncoupledPatchExpr,
            CoupledPatchExpr,
            value_type
        >
        (
            elems.dimensions(),
            elems.oriented()
        ),
        elems_(elems)
    {}

    //- Move construct
    GeometricFieldRefWrap(this_type&& elems)
    :
        GeometricFieldExpression
        <
            GeometricFieldRefWrap<GeoField>,
            IntExpr,
            UncoupledPatchExpr,
            CoupledPatchExpr,
            value_type
        >
        (
            elems.dimensions(),
            elems.oriented()
        ),
        elems_(elems)
    {}

    // Construct from ListExpression, forcing its evaluation.
    template<typename E>
    GeometricFieldRefWrap
    (
        this_type& elems,
        const GeometricFieldExpression
        <
            E,
            typename E::IntExpr,
            typename E::UncoupledPatchExpr,
            typename E::CoupledPatchExpr,
            typename E::value_type
        >& expr
    )
    :
        GeometricFieldExpression
        <
            GeometricFieldRefWrap<GeoField>,
            IntExpr,
            UncoupledPatchExpr,
            CoupledPatchExpr,
            value_type
        >
        (
            elems.dimensions(),
            elems.oriented()
        ),
        elems_(elems)
    {
        expr.evaluate(elems_);
    }

    //- Assignment
    template<typename E>
    void operator=
    (
        const GeometricFieldExpression
        <
            E,
            typename E::IntExpr,
            typename E::UncoupledPatchExpr,
            typename E::CoupledPatchExpr,
            typename E::value_type
        >& expr
    )
    {
        expr.evaluate(elems_);
    }

    //- Evaluate and return as GeoField. Rename to evaluate to make it clear
    //- it takes time? Or leave as indexing for convenience?
    template<typename E>
    this_type& evaluate
    (
        const GeometricFieldExpression
        <
            E,
            typename E::IntExpr,
            typename E::UncoupledPatchExpr,
            typename E::CoupledPatchExpr,
            typename E::value_type
        >& expr
    )
    {
        return expr.evaluate(elems_);
    }

    value_type operator[](const label i) const
    {
        return elems_[i];
    }

    value_type& operator[](const label  i)
    {
        return elems_[i];
    }

    auto size() const noexcept
    {
        return elems_.size();
    }

    IntExpr internalField()
    {
        auto& fld = elems_.internalFieldRef();
        return IntExpr(fld.size(), fld);
    }

    IntExpr internalField() const
    {
        auto& fld = elems_.internalField();
        return IntExpr(fld.size(), fld);
    }

    UncoupledPatchExpr patchField(const label i)
    {
        auto& fld = elems_.boundaryFieldRef()[i];
        if constexpr
        (
            std::is_same_v<pointMesh, typename GeoField::Mesh>
        )
        {
            const auto* vp = isA<Field<value_type>>(fld);
            if (vp)
            {
                auto& v = *vp;
                return UncoupledPatchExpr
                (
                    v.size(),
                    const_cast<Field<value_type>&>(v)
                );
            }
            else
            {
                return UncoupledPatchExpr(0, List<value_type>::null());
            }
        }
        else
        {
            return UncoupledPatchExpr(fld.size(), fld);
        }
    }

    UncoupledPatchExpr patchField(const label i) const
    {
        const auto& fld = elems_.boundaryField()[i];
        if constexpr
        (
            std::is_same_v<pointMesh, typename GeoField::Mesh>
        )
        {
            const auto* vp = isA<Field<value_type>>(fld);
            if (vp)
            {
                return UncoupledPatchExpr(vp->size(), *vp);
            }
            else
            {
                return UncoupledPatchExpr(0, List<value_type>::null());
            }
        }
        else
        {
            return UncoupledPatchExpr(fld.size(), fld);
        }
    }

    CoupledPatchExpr coupledPatchField(const label i)
    {
        auto& fld = elems_.boundaryFieldRef()[i];
        if constexpr
        (
            std::is_same_v<pointMesh, typename GeoField::Mesh>
        )
        {
            const auto* vp = isA<Field<value_type>>(fld);
            if (vp)
            {
                auto& v = *vp;
                return CoupledPatchExpr
                (
                    //v.size(),
                    const_cast<Field<value_type>&>(v)
                );
            }
            else
            {
                return CoupledPatchExpr(List<value_type>::null());
            }
        }
        else
        {
            return CoupledPatchExpr(fld);
        }
    }

    CoupledPatchExpr coupledPatchField(const label i) const
    {
        const auto& fld = elems_.boundaryField()[i];
        if constexpr
        (
            std::is_same_v<pointMesh, typename GeoField::Mesh>
        )
        {
            const auto* vp = isA<Field<value_type>>(fld);
            if (vp)
            {
                //return CoupledPatchExpr(vp->size(), *vp);
                return CoupledPatchExpr(*vp);
            }
            else
            {
                //return CoupledPatchExpr(0, List<value_type>::null());
                return CoupledPatchExpr(List<value_type>::null());
            }
        }
        else
        {
            //return UncoupledPatchExpr(fld.size(), fld);
            return CoupledPatchExpr(fld);
        }
    }

    template<class AccessOp>
    auto access(const AccessOp& cop, const label i) const
    {
        return cop(elems_, i);
    }
};


template<class GeoField>
class GeometricFieldConstRefWrap
:
    public GeometricFieldExpression
    <
        GeometricFieldConstRefWrap<GeoField>,
        ListConstRefWrap<typename GeoField::value_type>,
        ListConstRefWrap<typename GeoField::value_type>,
        ListConstTmpWrap<Field<typename GeoField::value_type>>,
        typename GeoField::value_type
    >
{
public:

    static constexpr bool is_leaf = false;  //true;

    //- The GeometricField type
    typedef GeoField this_type;

    //- Type to return for internal field
    typedef typename GeoField::value_type value_type;

    //- Type to return for patchField
    typedef ListConstRefWrap<value_type> IntExpr;
    typedef ListConstRefWrap<value_type> UncoupledPatchExpr;
    typedef ListConstTmpWrap<Field<value_type>> CoupledPatchExpr;


private:

    const this_type& elems_;


public:

    // Construct from components
    GeometricFieldConstRefWrap(const this_type& elems)
    :
        GeometricFieldExpression
        <
            GeometricFieldConstRefWrap<GeoField>,
            IntExpr,
            UncoupledPatchExpr,
            CoupledPatchExpr,
            value_type
        >
        (
            elems.dimensions(),
            elems.oriented()
        ),
        elems_(elems)
    {}

    // return the underlying data
    const auto& data() const
    {
        return elems_;
    }

    value_type operator[](const label i) const
    {
        return elems_[i];
    }

    auto size() const noexcept
    {
        return elems_.size();
    }

    IntExpr internalField() const
    {
        return elems_.internalField();
    }

    UncoupledPatchExpr patchField(const label i) const
    {
        //return elems_.boundaryField()[i];
        const auto& fld = elems_.boundaryField()[i];
        if constexpr
        (
            std::is_same_v<pointMesh, typename GeoField::Mesh>
        )
        {
            const auto* vp = isA<Field<value_type>>(fld);
            if (vp)
            {
                return UncoupledPatchExpr(*vp);
            }
            else
            {
                return UncoupledPatchExpr(List<value_type>::null());
            }
        }
        else
        {
            return UncoupledPatchExpr(fld);
        }
    }

    CoupledPatchExpr coupledPatchField(const label i) const
    {
        const auto& fld = elems_.boundaryField()[i];
        if constexpr
        (
            std::is_same_v<pointMesh, typename GeoField::Mesh>
        )
        {
            const auto* vp = isA<Field<value_type>>(fld);
            if (vp)
            {
                return CoupledPatchExpr(*vp);
            }
            else
            {
                return CoupledPatchExpr(tmp<Field<value_type>>(nullptr));
            }
        }
        else
        {
            return CoupledPatchExpr(fld);
        }
    }

    template<class AccessOp>
    auto access(const AccessOp& cop, const label i) const
    {
        return cop(elems_, i);
    }
};


template<class GeoField>
class GeometricFieldConstTmpWrap
:
    public GeometricFieldExpression
    <
        GeometricFieldConstTmpWrap<GeoField>,
        ListConstRefWrap<typename GeoField::value_type>,
        ListConstRefWrap<typename GeoField::value_type>,
        ListConstTmpWrap<Field<typename GeoField::value_type>>,
        typename GeoField::value_type
    >
{
public:

    //- Have expressions use copy to maintain tmp refCount
    static constexpr bool is_leaf = false;

    //- The GeometricField type
    typedef GeoField this_type;

    //- Type to return for internal field
    typedef typename GeoField::value_type value_type;

    //- Type to return for patchField
    typedef ListConstRefWrap<value_type> IntExpr;
    typedef ListConstRefWrap<value_type> UncoupledPatchExpr;
    typedef ListConstTmpWrap<Field<value_type>> CoupledPatchExpr;


private:

    const tmp<this_type> elems_;


public:

    // Construct from components
    GeometricFieldConstTmpWrap(const tmp<this_type>& elems)
    :
        GeometricFieldExpression
        <
            GeometricFieldConstTmpWrap<GeoField>,
            IntExpr,
            UncoupledPatchExpr,
            CoupledPatchExpr,
            value_type
        >
        (
            elems().dimensions(),
            elems().oriented()
        ),
        elems_(elems, true)
    {
        //Pout<< "GeometricFieldConstTmpWrap :" << nl
        //    << "    this_type  : "
        //    << boost::core::demangle(typeid(this_type).name())
        //    << nl
        //    << "    PatchExpr  : "
        //    << boost::core::demangle(typeid(IntExpr).name()) << nl
        //    << endl;
    }

    GeometricFieldConstTmpWrap
    (
        const GeometricFieldConstTmpWrap<GeoField>& w
    )
    :
        GeometricFieldExpression
        <
            GeometricFieldConstTmpWrap<GeoField>,
            IntExpr,
            UncoupledPatchExpr,
            CoupledPatchExpr,
            value_type
        >
        (
            w.dimensions(),
            w.oriented()
        ),
        //elems_(std::move(w.elems_))
        elems_(w.elems_)
    {
        //Pout<< "GeometricFieldConstTmpWrap copy :" << nl
        //    << "    this_type  : "
        //    << boost::core::demangle(typeid(this_type).name())
        //    << nl
        //    << "    PatchExpr  : "
        //    << boost::core::demangle(typeid(IntExpr).name()) << nl
        //    << endl;
    }

    ~GeometricFieldConstTmpWrap()
    {
        //Pout<< "~GeometricFieldConstTmpWrap destructor :" << nl
        //    << "    this_type  : "
        //    << boost::core::demangle(typeid(this_type).name())
        //    << nl
        //    << endl;
    }

    // return the underlying data
    const auto& data() const
    {
        return elems_();
    }

    value_type operator[](const label i) const
    {
        // Redirection for every access. Slow?
        return elems_()[i];
    }

    auto size() const noexcept
    {
        return elems_().size();
    }

    IntExpr internalField() const
    {
        return elems_().internalField();
    }

    UncoupledPatchExpr patchField(const label i) const
    {
        const auto& fld = elems_().boundaryField()[i];
        if constexpr
        (
            std::is_same_v<pointMesh, typename GeoField::Mesh>
        )
        {
            const auto* vp = isA<Field<value_type>>(fld);
            if (vp)
            {
                return UncoupledPatchExpr(*vp);
            }
            else
            {
                return UncoupledPatchExpr(List<value_type>::null());
            }
        }
        else
        {
            return UncoupledPatchExpr(fld);
        }

    }

    CoupledPatchExpr coupledPatchField(const label i) const
    {
        const auto& fld = elems_().boundaryField()[i];
        if constexpr
        (
            std::is_same_v<pointMesh, typename GeoField::Mesh>
        )
        {
            const auto* vp = isA<Field<value_type>>(fld);
            if (vp)
            {
                return CoupledPatchExpr(*vp);
            }
            else
            {
                return CoupledPatchExpr(List<value_type>::null());
            }
        }
        else
        {
            return CoupledPatchExpr(fld);
        }
    }

    template<class AccessOp>
    auto access(const AccessOp& cop, const label i) const
    {
        return cop(elems_(), i);
    }
};


// Expressions on GeometricFields
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Do '+' ourselves - just to show. Others done with macros below.
template<typename E1, typename E2>
class GF_add
:
    public GeometricFieldExpression
    <
        GF_add<E1, E2>,
        List_add
        <
            typename E1::IntExpr,
            typename E2::IntExpr
        >,
        List_add
        <
            typename E1::UncoupledPatchExpr,
            typename E2::UncoupledPatchExpr
        >,
        List_add
        <
            typename E1::CoupledPatchExpr,
            typename E2::CoupledPatchExpr
        >,
        typename E1::value_type
    >
{
    // cref if leaf, copy otherwise
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type u_;
    typename std::conditional<E2::is_leaf, const E2&, const E2>::type v_;

public:
    static constexpr bool is_leaf = false;

    //- Type to return for internal field
    typedef typename E1::value_type value_type;

    //- Type to return for patchField
    typedef List_add
    <
        typename E1::IntExpr,
        typename E2::IntExpr
    > IntExpr;
    typedef List_add
    <
        typename E1::UncoupledPatchExpr,
        typename E2::UncoupledPatchExpr
    > UncoupledPatchExpr;
    typedef List_add
    <
        typename E1::CoupledPatchExpr,
        typename E2::CoupledPatchExpr
    > CoupledPatchExpr;

    GF_add(const E1& u, const E2& v)
    :
        GeometricFieldExpression
        <
            GF_add<E1, E2>,
            IntExpr,
            UncoupledPatchExpr,
            CoupledPatchExpr,
            value_type
        >
        (
            u.dimensions() + v.dimensions(),
            u.oriented()    // + v.oriented()
        ),
        u_(u),
        v_(v)
    {
        assert(u.size() == v.size());
        static_assert
        (
            std::is_same_v
            <
                typename E1::value_type,
                typename E2::value_type
            > == true
        );
    }
    auto operator[](const label i) const
    {
        return u_[i] + v_[i];
    }
    auto size() const noexcept { return u_.size(); }

    IntExpr internalField() const
    {
        return IntExpr(u_.internalField(), v_.internalField());
    }

    UncoupledPatchExpr patchField(const label i) const
    {
        return UncoupledPatchExpr(u_.patchField(i), v_.patchField(i));
    }

    CoupledPatchExpr coupledPatchField(const label i) const
    {
        return CoupledPatchExpr
        (
            u_.coupledPatchField(i),
            v_.coupledPatchField(i)
        );
    }

    template<class AccessOp>
    auto access(const AccessOp& cop, const label i) const
    {
        return UncoupledPatchExpr(u_.access(cop, i), v_.access(cop, i));
    }
};
template<typename E1, typename E2>
GF_add<E1, E2>
operator+
(
    GeometricFieldExpression
    <
        E1,
        typename E1::IntExpr,
        typename E1::UncoupledPatchExpr,
        typename E1::CoupledPatchExpr,
        typename E1::value_type
    > const& u,
    GeometricFieldExpression
    <
        E2,
        typename E2::IntExpr,
        typename E2::UncoupledPatchExpr,
        typename E2::CoupledPatchExpr,
        typename E2::value_type
    > const& v
)
{
    return GF_add<E1, E2>
    (
        static_cast<const E1&>(u),
        static_cast<const E2&>(v)
    );
}


#define EXPRESSION_GF_FUNCTION1_DIMLESS(Func, BaseFunc, WrapType, OpFunc)      \
template<typename E1>                                                          \
class OpFunc                                                                   \
:                                                                              \
    public GeometricFieldExpression                                            \
    <                                                                          \
        OpFunc<E1>,                                                            \
        WrapType<typename E1::IntExpr>,                                        \
        WrapType<typename E1::UncoupledPatchExpr>,                             \
        WrapType<typename E1::CoupledPatchExpr>,                               \
        typename E1::value_type                                                \
    >                                                                          \
{                                                                              \
    /* cref if leaf, copy otherwise */                                         \
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type u_;      \
                                                                               \
public:                                                                        \
    static constexpr bool is_leaf = false;                                     \
                                                                               \
    /* Type to return for internal field */                                    \
    typedef typename E1::value_type value_type;                                \
                                                                               \
    /* Type to return for patchField */                                        \
    typedef WrapType<typename E1::IntExpr> IntExpr;                            \
    typedef WrapType<typename E1::UncoupledPatchExpr> UncoupledPatchExpr;      \
    typedef WrapType<typename E1::CoupledPatchExpr> CoupledPatchExpr;          \
                                                                               \
    OpFunc(const E1& u)                                                        \
    :                                                                          \
        GeometricFieldExpression                                               \
        <                                                                      \
            OpFunc<E1>,                                                        \
            IntExpr,                                                           \
            UncoupledPatchExpr,                                                \
            CoupledPatchExpr,                                                  \
            value_type                                                         \
        >                                                                      \
        (                                                                      \
            u.dimensions(),                                                    \
            u.oriented()                                                       \
        ),                                                                     \
        u_(u)                                                                  \
    {}                                                                         \
    auto operator[](const label i) const                                       \
    {                                                                          \
        return BaseFunc(u_[i]);                                                \
    }                                                                          \
    auto size() const noexcept { return u_.size(); }                           \
                                                                               \
    IntExpr internalField() const                                              \
    {                                                                          \
        return IntExpr(u_.internalField());                                    \
    }                                                                          \
                                                                               \
    UncoupledPatchExpr patchField(const label i) const                         \
    {                                                                          \
        return UncoupledPatchExpr(u_.patchField(i));                           \
    }                                                                          \
                                                                               \
    CoupledPatchExpr coupledPatchField(const label i) const                    \
    {                                                                          \
        return CoupledPatchExpr(u_.coupledPatchField(i));                      \
    }                                                                          \
                                                                               \
    template<class AccessOpOp>                                                 \
    auto access(const AccessOpOp& cop, const label i) const                    \
    {                                                                          \
        return UncoupledPatchExpr(u_.access(cop, i));                          \
    }                                                                          \
};                                                                             \
template<typename E1>                                                          \
OpFunc<E1> Func                                                                \
(                                                                              \
    GeometricFieldExpression                                                   \
    <                                                                          \
        E1,                                                                    \
        typename E1::IntExpr,                                                  \
        typename E1::UncoupledPatchExpr,                                       \
        typename E1::CoupledPatchExpr,                                         \
        typename E1::value_type                                                \
    > const& u                                                                 \
)                                                                              \
{                                                                              \
    return OpFunc<E1>(static_cast<const E1&>(u));                              \
}


#define EXPRESSION_GF_FUNCTION1(Func, BaseFunc, WrapType, OpFunc)              \
template<typename E1>                                                          \
class OpFunc                                                                   \
:                                                                              \
    public GeometricFieldExpression                                            \
    <                                                                          \
        OpFunc<E1>,                                                            \
        WrapType<typename E1::IntExpr>,                                        \
        WrapType<typename E1::UncoupledPatchExpr>,                             \
        WrapType<typename E1::CoupledPatchExpr>,                               \
        typename E1::value_type                                                \
    >                                                                          \
{                                                                              \
    /* cref if leaf, copy otherwise */                                         \
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type u_;      \
                                                                               \
public:                                                                        \
    static constexpr bool is_leaf = false;                                     \
                                                                               \
    /* Type to return for internal field */                                    \
    typedef typename E1::value_type value_type;                                \
                                                                               \
    /* Type to return for patchField */                                        \
    typedef WrapType<typename E1::IntExpr> IntExpr;                            \
    typedef WrapType<typename E1::UncoupledPatchExpr> UncoupledPatchExpr;      \
    typedef WrapType<typename E1::CoupledPatchExpr> CoupledPatchExpr;          \
                                                                               \
    OpFunc(const E1& u)                                                        \
    :                                                                          \
        GeometricFieldExpression                                               \
        <                                                                      \
            OpFunc<E1>,                                                        \
            IntExpr,                                                           \
            UncoupledPatchExpr,                                                \
            CoupledPatchExpr,                                                  \
            value_type                                                         \
        >                                                                      \
        (                                                                      \
            BaseFunc(u.dimensions()),                                          \
            BaseFunc(u.oriented())                                             \
        ),                                                                     \
        u_(u)                                                                  \
    {}                                                                         \
    auto operator[](const label i) const                                       \
    {                                                                          \
        return BaseFunc(u_[i]);                                                \
    }                                                                          \
    auto size() const noexcept { return u_.size(); }                           \
                                                                               \
    IntExpr internalField() const                                              \
    {                                                                          \
        return IntExpr(u_.internalField());                                    \
    }                                                                          \
                                                                               \
    UncoupledPatchExpr patchField(const label i) const                         \
    {                                                                          \
        return UncoupledPatchExpr(u_.patchField(i));                           \
    }                                                                          \
                                                                               \
    CoupledPatchExpr coupledPatchField(const label i) const                    \
    {                                                                          \
        return CoupledPatchExpr(u_.coupledPatchField(i));                      \
    }                                                                          \
                                                                               \
    template<class AccessOpOp>                                                 \
    auto access(const AccessOpOp& cop, const label i) const                    \
    {                                                                          \
        return UncoupledPatchExpr(u_.access(cop, i));                          \
    }                                                                          \
};                                                                             \
template<typename E1>                                                          \
OpFunc<E1> Func                                                                \
(                                                                              \
    GeometricFieldExpression                                                   \
    <                                                                          \
        E1,                                                                    \
        typename E1::IntExpr,                                                  \
        typename E1::UncoupledPatchExpr,                                       \
        typename E1::CoupledPatchExpr,                                         \
        typename E1::value_type                                                \
    > const& u                                                                 \
)                                                                              \
{                                                                              \
    return OpFunc<E1>(static_cast<const E1&>(u));                              \
}


//// Do negation separately - problem with macro expansion. Tbd.
//template<typename E1>
//class GeometricFieldNegate
//:
//    public GeometricFieldExpression
//    <
//        GeometricFieldNegate<E1>,
//        List_negate<typename E1::IntExpr>,
//        List_negate<typename E1::UncoupledPatchExpr>,
//        List_negate<typename E1::CoupledPatchExpr>,
//        typename E1::value_type
//    >
//{
//    /* cref if leaf, copy otherwise */ 
//    typename std::conditional<E1::is_leaf, const E1&, const E1>::type u_;
//
//public:
//    static constexpr bool is_leaf = false;
//
//    /* Type to return for internal field */
//    typedef typename E1::value_type value_type;
//
//    /* Type to return for patchField */
//    typedef List_negate<typename E1::IntExpr> IntExpr;
//    typedef List_negate<typename E1::UncoupledPatchExpr> UncoupledPatchExpr;
//    typedef List_negate<typename E1::CoupledPatchExpr> CoupledPatchExpr;
//
//    GeometricFieldNegate(const E1& u)
//    :
//        GeometricFieldExpression
//        <
//            GeometricFieldNegate<E1>,
//            IntExpr,
//            UncoupledPatchExpr,
//            CoupledPatchExpr,
//            value_type
//        >
//        (
//            -u.dimensions(),
//            -u.oriented()
//        ),
//        u_(u)
//    {}
//    auto operator[](const label i) const
//    {
//        return -u_[i];
//    }
//    label size() const { return u_.size(); }
//
//    IntExpr internalField() const
//    {
//        return PatchExpr(u_.internalField());
//    }
//
//    UncoupledPatchExpr patchField(const label i) const
//    {
//        return UncoupledPatchExpr(u_.patchField(i));
//    }
//};
//template<typename E1>
//GeometricFieldNegate<E1> operator-
//(
//    GeometricFieldExpression
//    <
//        E1,
//        typename E1::IntExpr,
//        typename E1::UncoupledPatchExpr,
//        typename E1::CoupledPatchExpr,
//        typename E1::value_type
//    > const& u
//)
//{
//    return GeometricFieldNegate<E1>(static_cast<const E1&>(u));
//}

#define EXPRESSION_GF_OPERATOR(Op, WrapType, OpFunc)                           \
template<typename E1, typename E2>                                             \
class OpFunc                                                                   \
:                                                                              \
    public GeometricFieldExpression                                            \
    <                                                                          \
        OpFunc<E1, E2>,                                                        \
        WrapType                                                               \
        <                                                                      \
            typename E1::IntExpr,                                              \
            typename E2::IntExpr                                               \
        >,                                                                     \
        WrapType                                                               \
        <                                                                      \
            typename E1::UncoupledPatchExpr,                                   \
            typename E2::UncoupledPatchExpr                                    \
        >,                                                                     \
        WrapType                                                               \
        <                                                                      \
            typename E1::CoupledPatchExpr,                                     \
            typename E2::CoupledPatchExpr                                      \
        >,                                                                     \
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
    /* Type to return for patchField */                                        \
    typedef WrapType                                                           \
    <                                                                          \
        typename E1::IntExpr,                                                  \
        typename E2::IntExpr                                                   \
    > IntExpr;                                                                 \
    typedef WrapType                                                           \
    <                                                                          \
        typename E1::UncoupledPatchExpr,                                       \
        typename E2::UncoupledPatchExpr                                        \
    > UncoupledPatchExpr;                                                      \
    typedef WrapType                                                           \
    <                                                                          \
        typename E1::CoupledPatchExpr,                                         \
        typename E2::CoupledPatchExpr                                          \
    > CoupledPatchExpr;                                                        \
                                                                               \
    OpFunc(const E1& u, const E2& v)                                           \
    :                                                                          \
        GeometricFieldExpression                                               \
        <                                                                      \
            OpFunc<E1, E2>,                                                    \
            IntExpr,                                                           \
            UncoupledPatchExpr,                                                \
            CoupledPatchExpr,                                                  \
            typename E1::value_type                                            \
        >                                                                      \
        (                                                                      \
            u.dimensions() Op v.dimensions(),                                  \
            u.oriented() Op v.oriented()                                       \
        ),                                                                     \
        u_(u),                                                                 \
        v_(v)                                                                  \
    {                                                                          \
      /*static_assert                                                          \
        (                                                                      \
            std::is_same_v                                                     \
            <                                                                  \
                typename E1::value_type,                                       \
                typename E2::value_type                                        \
            > == true                                                          \
        );*/                                                                   \
        assert(u.size() == v.size());                                          \
    }                                                                          \
    auto operator[](const label i) const                                       \
    {                                                                          \
        return u_[i] Op v_[i];                                                 \
    }                                                                          \
    auto size() const noexcept { return u_.size(); }                           \
                                                                               \
    IntExpr internalField() const                                              \
    {                                                                          \
        return IntExpr(u_.internalField(), v_.internalField());                \
    }                                                                          \
                                                                               \
    UncoupledPatchExpr patchField(const label i) const                         \
    {                                                                          \
        return UncoupledPatchExpr(u_.patchField(i), v_.patchField(i));         \
    }                                                                          \
                                                                               \
    CoupledPatchExpr coupledPatchField(const label i) const                    \
    {                                                                          \
        return CoupledPatchExpr                                                \
        (                                                                      \
            u_.coupledPatchField(i),                                           \
            v_.coupledPatchField(i)                                            \
        );                                                                     \
    }                                                                          \
                                                                               \
    template<class AccessOp>                                                   \
    auto access(const AccessOp& cop, const label i) const                      \
    {                                                                          \
        return UncoupledPatchExpr                                              \
        (                                                                      \
            u_.access(cop, i),                                                 \
            v_.access(cop, i)                                                  \
        );                                                                     \
    }                                                                          \
};                                                                             \
template<typename E1, typename E2>                                             \
OpFunc<E1, E2>                                                                 \
operator Op                                                                    \
(                                                                              \
    GeometricFieldExpression                                                   \
    <                                                                          \
        E1,                                                                    \
        typename E1::IntExpr,                                                  \
        typename E1::UncoupledPatchExpr,                                       \
        typename E1::CoupledPatchExpr,                                         \
        typename E1::value_type                                                \
    > const& u,                                                                \
    GeometricFieldExpression                                                   \
    <                                                                          \
        E2,                                                                    \
        typename E2::IntExpr,                                                  \
        typename E2::UncoupledPatchExpr,                                       \
        typename E2::CoupledPatchExpr,                                         \
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

#define EXPRESSION_GF_FUNCTION2(Func, BaseFunc, WrapType, OpFunc)              \
template<typename E1, typename E2>                                             \
class OpFunc                                                                   \
:                                                                              \
    public GeometricFieldExpression                                            \
    <                                                                          \
        OpFunc<E1, E2>,                                                        \
        WrapType                                                               \
        <                                                                      \
            typename E1::IntExpr,                                              \
            typename E2::IntExpr                                               \
        >,                                                                     \
        WrapType                                                               \
        <                                                                      \
            typename E1::UncoupledPatchExpr,                                   \
            typename E2::UncoupledPatchExpr                                    \
        >,                                                                     \
        WrapType                                                               \
        <                                                                      \
            typename E1::CoupledPatchExpr,                                     \
            typename E2::CoupledPatchExpr                                      \
        >,                                                                     \
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
    /* Type to return for patchField */                                        \
    typedef WrapType                                                           \
    <                                                                          \
        typename E1::IntExpr,                                                  \
        typename E2::IntExpr                                                   \
    > IntExpr;                                                                 \
    typedef WrapType                                                           \
    <                                                                          \
        typename E1::UncoupledPatchExpr,                                       \
        typename E2::UncoupledPatchExpr                                        \
    > UncoupledPatchExpr;                                                      \
    typedef WrapType                                                           \
    <                                                                          \
        typename E1::CoupledPatchExpr,                                         \
        typename E2::CoupledPatchExpr                                          \
    > CoupledPatchExpr;                                                        \
                                                                               \
    OpFunc(const E1& u, const E2& v)                                           \
    :                                                                          \
        GeometricFieldExpression                                               \
        <                                                                      \
            OpFunc<E1, E2>,                                                    \
            IntExpr,                                                           \
            UncoupledPatchExpr,                                                \
            CoupledPatchExpr,                                                  \
            value_type                                                         \
        >                                                                      \
        (                                                                      \
            BaseFunc(u.dimensions(), v.dimensions()),                          \
            BaseFunc(u.oriented(), v.oriented())                               \
        ),                                                                     \
        u_(u),                                                                 \
        v_(v)                                                                  \
    {                                                                          \
        /*static_assert                                                        \
        (                                                                      \
            std::is_same_v                                                     \
            <                                                                  \
                typename E1::value_type,                                       \
                typename E2::value_type                                        \
            > == true                                                          \
        );*/                                                                   \
    }                                                                          \
    auto operator[](const label i) const                                       \
    {                                                                          \
        return BaseFunc(u_[i], v_[i]);                                         \
    }                                                                          \
    auto size() const noexcept {return u_.size(); }                            \
                                                                               \
    IntExpr internalField() const                                              \
    {                                                                          \
        return IntExpr(u_.internalField(), v_.internalField());                \
    }                                                                          \
                                                                               \
    UncoupledPatchExpr patchField(const label i) const                         \
    {                                                                          \
        return UncoupledPatchExpr(u_.patchField(i), v_.patchField(i));         \
    }                                                                          \
                                                                               \
    CoupledPatchExpr coupledPatchField(const label i) const                    \
    {                                                                          \
        return CoupledPatchExpr                                                \
        (                                                                      \
            u_.coupledPatchField(i),                                           \
            v_.coupledPatchField(i)                                            \
        );                                                                     \
    }                                                                          \
                                                                               \
    template<class AccessOp>                                                   \
    auto access(const AccessOp& cop, const label i) const                      \
    {                                                                          \
        return UncoupledPatchExpr                                              \
        (                                                                      \
            u_.access(cop, i),                                                 \
            v_.access(cop, i)                                                  \
        );                                                                     \
    }                                                                          \
};                                                                             \
template<typename E1, typename E2>                                             \
OpFunc<E1, E2> Func                                                            \
(                                                                              \
    GeometricFieldExpression                                                   \
    <                                                                          \
        E1,                                                                    \
        typename E1::IntExpr,                                                  \
        typename E1::UncoupledPatchExpr,                                       \
        typename E1::CoupledPatchExpr,                                         \
        typename E1::value_type                                                \
    > const& u,                                                                \
    GeometricFieldExpression                                                   \
    <                                                                          \
        E2,                                                                    \
        typename E2::IntExpr,                                                  \
        typename E2::UncoupledPatchExpr,                                       \
        typename E2::CoupledPatchExpr,                                         \
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


//- GeometricField expressions.

// macro arguments:
//  - 'sin' : name of new function (in Expression namespace)
//  - '::sin' : per-element function to call
//  - 'List_sin' : helper class to handle patchFields
//  - 'GF_sin' : generated expression helper class
EXPRESSION_GF_FUNCTION1_DIMLESS(sin, ::sin, List_sin, GF_sin)
EXPRESSION_GF_FUNCTION1_DIMLESS(cos, ::cos, List_cos, GF_cos)
EXPRESSION_GF_FUNCTION1_DIMLESS(tan, ::tan, List_tan, GF_tan)
EXPRESSION_GF_FUNCTION1_DIMLESS(sinh, ::sinh, List_sinh, GF_sinh)
EXPRESSION_GF_FUNCTION1_DIMLESS(cosh, ::cosh, List_cosh, GF_cosh)
EXPRESSION_GF_FUNCTION1_DIMLESS(tanh, ::tanh, List_tanh, GF_tanh)

#undef EXPRESSION_GF_FUNCTION1_DIMLESS

EXPRESSION_GF_FUNCTION1(sqr, Foam::sqr, List_sqr, GF_sqr)
EXPRESSION_GF_FUNCTION1(sqrt, Foam::sqrt, List_sqrt, GF_sqrt)
EXPRESSION_GF_FUNCTION1(magSqr, Foam::magSqr, List_magSqr, GF_magSqr)

EXPRESSION_GF_FUNCTION1(mag, Foam::mag, List_mag, GF_mag)
EXPRESSION_GF_FUNCTION1(symm, Foam::symm, List_symm, GF_symm)
EXPRESSION_GF_FUNCTION1(pow2, Foam::pow2, List_pow2, GF_pow2)
EXPRESSION_GF_FUNCTION1(pow3, Foam::pow3, List_pow3, GF_pow3)
EXPRESSION_GF_FUNCTION1(pow4, Foam::pow4, List_pow4, GF_pow4)

#undef EXPRESSION_GF_FUNCTION1

EXPRESSION_GF_FUNCTION2(min, Foam::min, List_min, GF_min)
EXPRESSION_GF_FUNCTION2(max, Foam::max, List_max, GF_max)

#undef EXPRESSION_GF_FUNCTION2

EXPRESSION_GF_OPERATOR(-, List_subtract, GF_subtract)
EXPRESSION_GF_OPERATOR(*, List_multiply, GF_multiply)
EXPRESSION_GF_OPERATOR(/, List_divide, GF_divide)
EXPRESSION_GF_OPERATOR(&, List_dot, GF_dot)

#undef EXPRESSION_GF_OPERATOR



// Expressions on constants
// ~~~~~~~~~~~~~~~~~~~~~~~~

template<class GeoField>
class UniformGeometricFieldWrap
:
    public GeometricFieldExpression
    <
        UniformGeometricFieldWrap<GeoField>,
        UniformListWrap<typename GeoField::value_type>,
        UniformListWrap<typename GeoField::value_type>,
        UniformListWrap<typename GeoField::value_type>,
        typename GeoField::value_type
    >
{
public:

    static constexpr bool is_leaf = false;  //true;

    //- Type to return for internal field
    typedef typename GeoField::value_type value_type;

    //- Type to return for patchField
    typedef UniformListWrap<value_type> IntExpr;
    typedef UniformListWrap<value_type> UncoupledPatchExpr;
    typedef UniformListWrap<value_type> CoupledPatchExpr;


private:

    const GeoField& elems_;

    const value_type val_;

public:

    // Construct from field (sizes only) and value
    UniformGeometricFieldWrap(const GeoField& elems, const value_type val)
    :
        GeometricFieldExpression
        <
            UniformGeometricFieldWrap<GeoField>,
            IntExpr,
            UncoupledPatchExpr,
            CoupledPatchExpr,
            value_type
        >
        (
            dimless,
            orientedType()
        ),
        elems_(elems),
        val_(val)
    {}

    // Construct from field (sizes only) and value, dimension
    UniformGeometricFieldWrap
    (
        const GeoField& elems,
        const dimensioned<value_type> val
    )
    :
        GeometricFieldExpression
        <
            UniformGeometricFieldWrap<GeoField>,
            IntExpr,
            UncoupledPatchExpr,
            CoupledPatchExpr,
            value_type
        >
        (
            val.dimensions(),
            orientedType()
        ),
        elems_(elems),
        val_(val.value())
    {}

    value_type operator[](const label i) const
    {
        return val_;
    }

    auto size() const noexcept
    {
        return elems_.size();
    }

    IntExpr internalField() const
    {
        return IntExpr(elems_.internalField().size(), val_);
    }

    UncoupledPatchExpr patchField(const label i) const
    {
        return UncoupledPatchExpr(elems_.boundaryField()[i].size(), val_);
    }

    CoupledPatchExpr coupledPatchField(const label i) const
    {
        return CoupledPatchExpr(elems_.boundaryField()[i].size(), val_);
    }

    template<class AccessOp>
    auto access(const AccessOp& cop, const label i) const
    {
        // Size+value only
        return patchField(i);
    }
};


//- Fully self-contained constant field wrapper. Not needed?
//template<class GeoField>
//class UniformGeometricFieldWrap2
//:
//    public GeometricFieldExpression
//    <
//        UniformGeometricFieldWrap2<GeoField>,
//        UniformListWrap<typename GeoField::value_type>,
//        typename GeoField::value_type
//    >
//{
//public:
//
//    static constexpr bool is_leaf = true;
//
//    //- Type to return for internal field
//    typedef typename GeoField::value_type value_type;
//
//    //- Type to return for patchField
//    typedef UniformListWrap<value_type> PatchExpr;
//
//
//private:
//
//    const value_type val_;
//
//    const label nInternal_;
//
//    labelList patchSizes_;
//
//    //const dimensionSet& dimensions_;
//
//public:
//
//    // Construct from field (sizes only) and value
//    UniformGeometricFieldWrap2(const GeoField& elems, const value_type val)
//    :
//        val_(val),
//        nInternal_(elems.size()),
//        patchSizes_(elems.boundaryField().size())
//        //dimensions_(elems.dimensions()),
//    {
//        for (const auto& pfld : elems.boundaryField())
//        {
//            patchSizes_[pfld.patch().index()] = pfld.size();
//        }
//    }
//
//    // Construct from field (sizes only) and value, dimension
//    UniformGeometricFieldWrap2
//    (
//        const GeoField& elems,
//        const dimensioned<value_type> val
//    )
//    :
//        val_(val.value()),
//        nInternal_(elems.size()),
//        patchSizes_(elems.boundaryField().size())
//        //dimensions_(val.dimensions()),
//    {
//        for (const auto& pfld : elems.boundaryField())
//        {
//            patchSizes_[pfld.patch().index()] = pfld.size();
//        }
//    }
//
//    value_type operator[](const label i) const
//    {
//        return val_;
//    }
//
//    auto size() const noexcept
//    {
//        return nInternal_;
//    }
//
//    PatchExpr internalField() const
//    {
//        return PatchExpr(nInternal_, val_);
//    }
//
//    PatchExpr patchField(const label i) const
//    {
//        return PatchExpr(patchSizes_[i], val_);
//    }
//};


/*
// Some fvm functions using expressions
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Supply a multiplier, either 1 or -1
template<class Expr, class Expr2>
void Su
(
    fvMatrix<typename Expr::value_type>& m,
    const Expr2& mult,
    const Expr& expression
)
{
    //- Wrap of List as an expression
    typedef ListConstRefWrap<typename Expr::value_type> expr;
    //- Evaluator of an expression
    typedef ListRefWrap<typename Expr::value_type> evaluator;

    // Wrap mesh volume
    const expr V(m.psi().mesh().V());

    // Add expression to source
    auto& s = m.source();
    evaluator(s, expr(s) + mult*V*expression);
}


// Always add
template<class Expr>
void Su
(
    fvMatrix<typename Expr::value_type>& m,
    const Expr& expression
)
{
    //- Wrap of List as an expression
    typedef ListConstRefWrap<typename Expr::value_type> expr;
    //- Evaluator of an expression
    typedef ListRefWrap<typename Expr::value_type> evaluator;

    // Wrap mesh volume
    const expr V(m.psi().mesh().V());

    // Add expression to source
    auto& s = m.source();
    evaluator(s, expr(s) + V*expression);
}


template<class Expr>
void rhs
(
    fvMatrix<typename Expr::value_type>& m,
    const Expr& expression
)
{
    Su(m, expression);
}


template<class Expr, class Expr2>
void Sp
(
    fvMatrix<typename Expr::value_type>& m,
    const Expr2& mult,
    const Expr& expression
)
{
    //- Wrap of List as an expression
    typedef ListConstRefWrap<typename Expr::value_type> expr;
    //- Evaluator of an expression
    typedef ListRefWrap<typename Expr::value_type> evaluator;

    // Wrap mesh volume
    const expr V(m.psi().mesh().V());

    // Add expression onto diag
    auto& d = m.diag();
    evaluator(d, expr(d) - mult*V*expression);
}


template<class Expr, class Expr2>
void SuSp
(
    fvMatrix<typename Expr::value_type>& m,
    const Expr2& mult,
    const Expr& expression
)
{
    //- Wrap of constant as a list expression
    typedef UniformListWrap<scalar> constant;
    //- Wrap of List as an expression
    typedef ListConstRefWrap<typename Expr::value_type> expr;
    //- Evaluator of an expression
    typedef ListRefWrap<typename Expr::value_type> evaluator;

    const constant constantZero(expression.size(), 0.0);

    // Wrap mesh volume
    const expr V(m.psi().mesh().V());

    // Linearise expression
    auto& diag = m.diag();
    auto& source = m.source();

    // diag() += domain*max(susp.field(), scalar(0));
    evaluator
    (
        diag,
        expr(diag) - mult*V*max(expression, constantZero)
    );
    // source() -= domain*min(susp.field(), scalar(0))*fld.primitiveField();
    evaluator
    (
        source,
        expr(source)
      + (
            mult
           *V
           *min(expression, constantZero)
           *expr(m.psi().internalField())
        )
    );
}
*/


} // End namespace Expression

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
