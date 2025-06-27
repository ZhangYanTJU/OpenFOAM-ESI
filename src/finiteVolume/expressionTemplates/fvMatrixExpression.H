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
    Foam::Expression

Description
    Expression templates for fvMatrix.

    TBD:
    - useImplicit_ etc.
    - lduMatrix::lowerCSR ?

SourceFiles
    fvMatrixExpression.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fvMatrixExpression_H
#define Foam_fvMatrixExpression_H

#include "ListExpression.H"
#include "GeometricFieldExpression.H"
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

// Expression of fvMatrix
// ~~~~~~~~~~~~~~~~~~~~~~
// Bit different:
// - immediately stores state (hasDiag, hasUpper etc.)
// - delays execution of building coefficients
// - big problem: individual components can be built differently so
//   different template argument ... Needs to be templated on all outputs
//   separately!
template
<
    typename E,
    typename DiagExpr,
    typename UpperExpr,
    typename LowerExpr,
    typename FaceFluxExpr,
    typename SourceExpr
>
class fvMatrixExpression
{
protected:

    const bool hasDiag_;
    const bool hasUpper_;
    const bool hasLower_;
    const bool hasFaceFluxCorrection_;
    dimensionSet dimensions_;

public:
    static constexpr bool is_leaf = false;

    // Constructors

        fvMatrixExpression
        (
            const bool hasDiag,
            const bool hasUpper,
            const bool hasLower,
            const bool hasFaceFluxCorrection,
            const dimensionSet& dimensions
        )
        :
            hasDiag_(hasDiag),
            hasUpper_(hasUpper),
            hasLower_(hasLower),
            hasFaceFluxCorrection_(hasFaceFluxCorrection),
            dimensions_(dimensions)
        {
            //Pout<< "fvMatrixExpression :" << nl
            //    << "    E        : "
            //    << boost::core::demangle(typeid(E).name())
            //    << nl
            //    << "    DiagExpr : "
            //    << boost::core::demangle(typeid(DiagExpr).name())
            //    << endl;
        }


    // lduMatrix part

        DiagExpr diag() const
        {
            return static_cast<const E&>(*this).diag();
        }

        DiagExpr diag()
        {
            return static_cast<E&>(*this).diag();
        }

        UpperExpr upper() const
        {
            return static_cast<const E&>(*this).upper();
        }

        UpperExpr upper()
        {
            return static_cast<E&>(*this).upper();
        }

        LowerExpr lower() const
        {
            return static_cast<const E&>(*this).lower();
        }

        LowerExpr lower()
        {
            return static_cast<E&>(*this).lower();
        }

        bool hasDiag() const
        {
             return hasDiag_;
        }
        bool hasUpper() const
        {
             return hasUpper_;
        }
        bool hasLower() const
        {
             return hasLower_;
        }
        bool diagonal() const
        {
             return (hasDiag() && !hasLower() && !hasUpper());
        }
        bool symmetric() const
        {
             return (hasDiag() && !hasLower() && hasUpper());
        }
        bool asymmetric() const
        {
             return (hasDiag() && hasLower() && hasUpper());
        }
        bool hasFaceFluxCorrection() const
        {
             return hasFaceFluxCorrection_;
        }


    // fvMatrix

        SourceExpr source()
        {
            return static_cast<E&>(*this).source();
        }

        SourceExpr source() const
        {
            return static_cast<const E&>(*this).source();
        }

        UpperExpr internalCoeffs(const label i)
        {
            return static_cast<E&>(*this).internalCoeffs(i);
        }

        UpperExpr internalCoeffs(const label i) const
        {
            return static_cast<const E&>(*this).internalCoeffs(i);
        }

        UpperExpr boundaryCoeffs(const label i)
        {
            return static_cast<E&>(*this).boundaryCoeffs(i);
        }

        UpperExpr boundaryCoeffs(const label i) const
        {
            return static_cast<const E&>(*this).boundaryCoeffs(i);
        }

        FaceFluxExpr faceFluxCorrection() const
        {
            return static_cast<const E&>(*this).faceFluxCorrection();
        }

        FaceFluxExpr faceFluxCorrection()
        {
            return static_cast<E&>(*this).faceFluxCorrection();
        }

        const dimensionSet& dimensions() const noexcept
        {
            return dimensions_;
        }


    // Other

        //- Helper to evaluate (=construct) an fvMatrix
        template<class Matrix>
        Matrix& evaluate(Matrix& m) const
        {
            if (hasDiag())
            {
                // Assign to diag. Also sets hasDiag.
                diag().evaluate(m.diag());
            }
            if (hasUpper())
            {
                upper().evaluate(m.upper());
            }
            if (hasLower())
            {
                lower().evaluate(m.lower());
            }
            if (hasFaceFluxCorrection())
            {
                faceFluxCorrection().evaluate(*m.faceFluxCorrectionPtr());
            }

            // Do boundary
            auto& intCoeffs = m.internalCoeffs();
            auto& bouCoeffs = m.boundaryCoeffs();
            const label n = intCoeffs.size();
            for (label i = 0; i < n; ++i)
            {
                if (intCoeffs.set(i) && intCoeffs[i].size())
                {
                    const auto intExpr = internalCoeffs(i);
                    intExpr.evaluate(intCoeffs[i]);
                }
                if (bouCoeffs.set(i) && bouCoeffs[i].size())
                {
                    const auto bouExpr = boundaryCoeffs(i);
                    bouExpr.evaluate(bouCoeffs[i]);
                }
            }

            const_cast<dimensionSet&>(m.dimensions()) = dimensions();

            return m;
        }
};


template<class Matrix>
class fvMatrixRefWrap
:
    public fvMatrixExpression
    <
        fvMatrixRefWrap<Matrix>,
        ListRefWrap<typename Matrix::psiFieldType::value_type>, // diag
        ListRefWrap<typename Matrix::psiFieldType::value_type>, // upper
        ListRefWrap<typename Matrix::psiFieldType::value_type>, // lower
        //ListRefWrap<typename Matrix::psiFieldType::value_type>, // faceflux
        GeometricFieldRefWrap<typename Matrix::faceFluxFieldType>,
        ListRefWrap<typename Matrix::psiFieldType::value_type>  // source
    >
{
public:

    static constexpr bool is_leaf = false;  //true;

    //- The fvMatrix type
    typedef Matrix this_type;

    //- Type to return for internal field
    typedef typename Matrix::psiFieldType::value_type value_type;

    //- Type to return for containers
    typedef ListRefWrap<value_type> DiagExpr;
    typedef ListRefWrap<value_type> UpperExpr;
    typedef ListRefWrap<value_type> LowerExpr;
    //typedef ListRefWrap<value_type> FaceFluxExpr;
    typedef GeometricFieldRefWrap
    <
        typename Matrix::faceFluxFieldType
    > FaceFluxExpr;
    typedef ListRefWrap<value_type> SourceExpr;


private:

    this_type& elems_;


public:

    //- Copy construct
    fvMatrixRefWrap(this_type& elems)
    :
        fvMatrixExpression
        <
            fvMatrixRefWrap<Matrix>,
            DiagExpr,
            UpperExpr,
            LowerExpr,
            FaceFluxExpr,
            SourceExpr
        >
        (
            elems.hasDiag(),
            elems.hasUpper(),
            elems.hasLower(),
            elems.hasFaceFluxCorrection(),
            elems.dimensions()
        ),
        elems_(elems)
    {}

    //- Move construct
    fvMatrixRefWrap(this_type&& elems)
    :
        fvMatrixExpression
        <
            fvMatrixRefWrap<Matrix>,
            DiagExpr,
            UpperExpr,
            LowerExpr,
            FaceFluxExpr,
            SourceExpr
        >
        (
            elems.hasDiag(),
            elems.hasUpper(),
            elems.hasLower(),
            elems.hasFaceFluxCorrection(),
            elems.dimensions()
        ),
        elems_(elems)
    {}

    // Construct from Matrix, forcing its evaluation.
    template<typename E>
    fvMatrixRefWrap
    (
        this_type& elems,
        const fvMatrixExpression
        <
            E,
            typename E::DiagExpr,
            typename E::UpperExpr,
            typename E::LowerExpr,
            typename E::FaceFluxExpr,
            typename E::SourceExpr
        >& expr
    )
    :
        fvMatrixExpression
        <
            fvMatrixRefWrap<Matrix>,
            DiagExpr,
            UpperExpr,
            LowerExpr,
            FaceFluxExpr,
            SourceExpr
        >
        (
            expr.hasDiag(),
            expr.hasUpper(),
            expr.hasLower(),
            elems.hasFaceFluxCorrection(),
            elems.dimensions()
        ),
        elems_(elems)
    {
        expr.evaluate(elems_);
    }

    //- Assignment
    template<typename E>
    void operator=
    (
        const fvMatrixExpression
        <
            E,
            typename E::DiagExpr,
            typename E::UpperExpr,
            typename E::LowerExpr,
            typename E::FaceFluxExpr,
            typename E::SourceExpr
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
        const fvMatrixExpression
        <
            E,
            typename E::DiagExpr,
            typename E::UpperExpr,
            typename E::LowerExpr,
            typename E::FaceFluxExpr,
            typename E::SourceExpr
        >& expr
    )
    {
        return expr.evaluate(elems_);
    }


    // lduMatrix part

        DiagExpr diag() const
        {
            return elems_.diag();
        }

        DiagExpr diag()
        {
            return elems_.diag();
        }

        UpperExpr upper() const
        {
            return elems_.upper();
        }

        UpperExpr upper()
        {
            return elems_.upper();
        }

        LowerExpr lower() const
        {
            return elems_.lower();
        }

        LowerExpr lower()
        {
            return elems_.lower();
        }


    // fvMatrix

        SourceExpr source() const
        {
            return elems_.source();
        }

        SourceExpr source()
        {
            return elems_.source();
        }

        UpperExpr internalCoeffs(const label i)
        {
            return
            (
                elems_.internalCoeffs().set(i)
              ? elems_.internalCoeffs()[i]
              : List<value_type>::null()
            );
        }

        UpperExpr internalCoeffs(const label i) const
        {
            return
            (
                elems_.internalCoeffs().set(i)
              ? elems_.internalCoeffs()[i]
              : List<value_type>::null()
            );
        }

        UpperExpr boundaryCoeffs(const label i)
        {
            return
            (
                elems_.boundaryCoeffs().set(i)
              ? elems_.boundaryCoeffs()[i]
              : List<value_type>::null()
            );
        }

        UpperExpr boundaryCoeffs(const label i) const
        {
            return
            (
                elems_.boundaryCoeffs().set(i)
              ? elems_.boundaryCoeffs()[i]
              : List<value_type>::null()
            );
        }

        FaceFluxExpr faceFluxCorrection() const
        {
            return *const_cast<this_type&>
            (
                elems_
            ).faceFluxCorrectionPtr().get();
        }

        FaceFluxExpr faceFluxCorrection()
        {
            return *elems_.faceFluxCorrectionPtr().get();
        }
};


template<class Matrix>
class fvMatrixConstRefWrap
:
    public fvMatrixExpression
    <
        fvMatrixConstRefWrap<Matrix>,
        ListConstRefWrap<typename Matrix::psiFieldType::value_type>, // diag
        ListConstRefWrap<typename Matrix::psiFieldType::value_type>, // upper
        ListConstRefWrap<typename Matrix::psiFieldType::value_type>, // lower
        //ListConstRefWrap<typename Matrix::psiFieldType::value_type>, // faceflux
        GeometricFieldConstRefWrap<typename Matrix::faceFluxFieldType>,
        ListConstRefWrap<typename Matrix::psiFieldType::value_type>  // source
    >
{
public:

    static constexpr bool is_leaf = false;  //true;

    //- The fvMatrix type
    typedef Matrix this_type;

    //- Type to return for internal field
    typedef typename Matrix::psiFieldType::value_type value_type;

    //- Type to return for containers
    typedef ListConstRefWrap<value_type> DiagExpr;
    typedef ListConstRefWrap<value_type> UpperExpr;
    typedef ListConstRefWrap<value_type> LowerExpr;
    //typedef ListConstRefWrap<value_type> FaceFluxExpr;
    typedef GeometricFieldConstRefWrap
    <
        typename Matrix::faceFluxFieldType
    > FaceFluxExpr;
    typedef ListConstRefWrap<value_type> SourceExpr;


private:

    const this_type& elems_;


public:

    // Construct from components
    fvMatrixConstRefWrap(const this_type& elems)
    :
        fvMatrixExpression
        <
            fvMatrixConstRefWrap<Matrix>,
            DiagExpr,
            UpperExpr,
            LowerExpr,
            FaceFluxExpr,
            SourceExpr
        >
        (
            elems.hasDiag(),
            elems.hasUpper(),
            elems.hasLower(),
            elems.hasFaceFluxCorrection(),
            elems.dimensions()
        ),
        elems_(elems)
    {}


    // lduMatrix part

        DiagExpr diag() const
        {
            return elems_.diag();
        }

        UpperExpr upper() const
        {
            return elems_.upper();
        }

        LowerExpr lower() const
        {
            return elems_.lower();
        }


    // fvMatrix

        SourceExpr source() const
        {
            return elems_.source();
        }

        UpperExpr internalCoeffs(const label i) const
        {
            return
            (
                elems_.internalCoeffs().set(i)
              ? elems_.internalCoeffs()[i]
              : List<value_type>::null()
            );
        }

        UpperExpr boundaryCoeffs(const label i) const
        {
            return
            (
                elems_.boundaryCoeffs().set(i)
              ? elems_.boundaryCoeffs()[i]
              : List<value_type>::null()
            );
        }

        FaceFluxExpr faceFluxCorrection() const
        {
            return *const_cast<this_type&>
            (
                elems_
            ).faceFluxCorrectionPtr().get();
        }
};


// Expressions on fvMatrices
// ~~~~~~~~~~~~~~~~~~~~~~~~~

// operator+
template<typename E1, typename E2>
class fvm_add
:
    public fvMatrixExpression
    <
        fvm_add<E1, E2>,
        List_add                    // DiagExpr
        <
            typename E1::DiagExpr,
            typename E2::DiagExpr
        >,
        List_add                    // UpperExpr
        <
            typename E1::UpperExpr,
            typename E2::UpperExpr
        >,
        List_add                    // LowerExpr
        <
            typename E1::LowerExpr,
            typename E2::LowerExpr
        >,
        GF_add                    // FaceFluxExpr
        <
            typename E1::FaceFluxExpr,
            typename E2::FaceFluxExpr
        >,
        List_add                    // SourceExpr
        <
            typename E1::SourceExpr,
            typename E2::SourceExpr
        >
    >
{
    // cref if leaf, copy otherwise
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type u_;
    typename std::conditional<E2::is_leaf, const E2&, const E2>::type v_;

public:
    static constexpr bool is_leaf = false;

    //- Type to return for internal field
    typedef typename E1::value_type value_type;

    //- Type to return for containers
    typedef List_add
    <
        typename E1::DiagExpr,
        typename E2::DiagExpr
    > DiagExpr;
    typedef List_add
    <
        typename E1::UpperExpr,
        typename E2::UpperExpr
    > UpperExpr;
    typedef List_add
    <
        typename E1::LowerExpr,
        typename E2::LowerExpr
    > LowerExpr;
    typedef GF_add
    <
        typename E1::FaceFluxExpr,
        typename E2::FaceFluxExpr
    > FaceFluxExpr;
    typedef List_add
    <
        typename E1::SourceExpr,
        typename E2::SourceExpr
    > SourceExpr;

    fvm_add(const E1& u, const E2& v)
    :
        fvMatrixExpression
        <
            fvm_add<E1, E2>,
            DiagExpr,
            UpperExpr,
            LowerExpr,
            FaceFluxExpr,
            SourceExpr
        >
        (
            u.hasDiag() || v.hasDiag(),
            u.hasUpper() || v.hasUpper(),
            u.hasLower() || v.hasLower(),
            u.hasFaceFluxCorrection() || v.hasFaceFluxCorrection(),
            u.dimensions()
        ),
        u_(u),
        v_(v)
    {
        //Pout<< "fvm_add :" << nl
        //    << "    u_ : " << boost::core::demangle(typeid(u_).name())
        //    << nl
        //    << "    v_ : " << boost::core::demangle(typeid(v_).name())
        //    << endl;
    }


    // lduMatrix part

        DiagExpr diag() const
        {
            return u_.diag() + v_.diag();
        }

        UpperExpr upper() const
        {
            return u_.upper() + v_.upper();
        }

        LowerExpr lower() const
        {
            return u_.lower() + v_.lower();
        }

    // fvMatrix

        SourceExpr source() const
        {
            return u_.source() + v_.source();
        }

        UpperExpr boundaryCoeffs(const label i) const
        {
            return u_.boundaryCoeffs(i) + v_.boundaryCoeffs(i);
        }

        UpperExpr internalCoeffs(const label i) const
        {
            return u_.internalCoeffs(i) + v_.internalCoeffs(i);
        }

        FaceFluxExpr faceFluxCorrection() const
        {
            return u_.faceFluxCorrection() + v_.faceFluxCorrection();
        }
};
template<typename E1, typename E2>
fvm_add<E1, E2>
operator+
(
    fvMatrixExpression
    <
        E1,
        typename E1::DiagExpr,
        typename E1::UpperExpr,
        typename E1::LowerExpr,
        typename E1::FaceFluxExpr,
        typename E1::SourceExpr
    > const& u,
    fvMatrixExpression
    <
        E2,
        typename E2::DiagExpr,
        typename E2::UpperExpr,
        typename E2::LowerExpr,
        typename E2::FaceFluxExpr,
        typename E2::SourceExpr
    > const& v
)
{
    return fvm_add<E1, E2>
    (
        static_cast<const E1&>(u),
        static_cast<const E2&>(v)
    );
}


// operator-, operator==
template<typename E1, typename E2>
class fvm_subtract
:
    public fvMatrixExpression
    <
        fvm_subtract<E1, E2>,
        List_subtract                   // diag
        <
            typename E1::DiagExpr,
            typename E2::DiagExpr
        >,
        List_subtract                   // upper
        <
            typename E1::UpperExpr,
            typename E2::UpperExpr
        >,
        List_subtract                   // lower
        <
            typename E1::LowerExpr,
            typename E2::LowerExpr
        >,
        GF_subtract                   // faceflux
        <
            typename E1::FaceFluxExpr,
            typename E2::FaceFluxExpr
        >,
        List_subtract                   // source
        <
            typename E1::SourceExpr,
            typename E2::SourceExpr
        >
    >
{
    // cref if leaf, copy otherwise
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type u_;
    typename std::conditional<E2::is_leaf, const E2&, const E2>::type v_;

public:
    static constexpr bool is_leaf = false;

    //- Type to return for internal field
    typedef typename E1::value_type value_type;

    //- Type to return for containers
    typedef List_subtract
    <
        typename E1::DiagExpr,
        typename E2::DiagExpr
    > DiagExpr;
    typedef List_subtract
    <
        typename E1::UpperExpr,
        typename E2::UpperExpr
    > UpperExpr;
    typedef List_subtract
    <
        typename E1::LowerExpr,
        typename E2::LowerExpr
    > LowerExpr;
    typedef GF_subtract
    <
        typename E1::FaceFluxExpr,
        typename E2::FaceFluxExpr
    > FaceFluxExpr;
    typedef List_subtract
    <
        typename E1::SourceExpr,
        typename E2::SourceExpr
    > SourceExpr;

    fvm_subtract(const E1& u, const E2& v)
    :
        fvMatrixExpression
        <
            fvm_subtract<E1, E2>,
            DiagExpr,
            UpperExpr,
            LowerExpr,
            FaceFluxExpr,
            SourceExpr
        >
        (
            u.hasDiag() || v.hasDiag(),
            u.hasUpper() || v.hasUpper(),
            u.hasLower() || v.hasLower(),
            u.hasFaceFluxCorrection() || v.hasFaceFluxCorrection(),
            u.dimensions()
        ),
        u_(u),
        v_(v)
    {}


    // lduMatrix part

        DiagExpr diag() const
        {
            return u_.diag() - v_.diag();
        }

        UpperExpr upper() const
        {
            return u_.upper() - v_.upper();
        }

        LowerExpr lower() const
        {
            return u_.lower() - v_.lower();
        }


    // fvMatrix

        SourceExpr source() const
        {
            return u_.source() - v_.source();
        }

        UpperExpr internalCoeffs(const label i) const
        {
            return u_.internalCoeffs(i) - v_.internalCoeffs(i);
        }

        UpperExpr boundaryCoeffs(const label i) const
        {
            return u_.boundaryCoeffs(i) - v_.boundaryCoeffs(i);
        }

        FaceFluxExpr faceFluxCorrection() const
        {
            return u_.faceFluxCorrection() - v_.faceFluxCorrection();
        }
};
template<typename E1, typename E2>
fvm_subtract<E1, E2>
operator-
(
    fvMatrixExpression
    <
        E1,
        typename E1::DiagExpr,
        typename E1::UpperExpr,
        typename E1::LowerExpr,
        typename E1::FaceFluxExpr,
        typename E1::SourceExpr
    > const& u,
    fvMatrixExpression
    <
        E2,
        typename E2::DiagExpr,
        typename E2::UpperExpr,
        typename E2::LowerExpr,
        typename E2::FaceFluxExpr,
        typename E2::SourceExpr
    > const& v
)
{
    return fvm_subtract<E1, E2>
    (
        static_cast<const E1&>(u),
        static_cast<const E2&>(v)
    );
}
template<typename E1, typename E2>
fvm_subtract<E1, E2>
operator==
(
    fvMatrixExpression
    <
        E1,
        typename E1::DiagExpr,
        typename E1::UpperExpr,
        typename E1::LowerExpr,
        typename E1::FaceFluxExpr,
        typename E1::SourceExpr
    > const& u,
    fvMatrixExpression
    <
        E2,
        typename E2::DiagExpr,
        typename E2::UpperExpr,
        typename E2::LowerExpr,
        typename E2::FaceFluxExpr,
        typename E2::SourceExpr
    > const& v
)
{
    return fvm_subtract<E1, E2>
    (
        static_cast<const E1&>(u),
        static_cast<const E2&>(v)
    );
}


// Negation : operator-
template<typename E1>
class fvm_negate
:
    public fvMatrixExpression
    <
        fvm_negate<E1>,
        List_negate<typename E1::DiagExpr>,
        List_negate<typename E1::UpperExpr>,
        List_negate<typename E1::LowerExpr>,
        List_negate<typename E1::FaceFluxExpr>,
        List_negate<typename E1::SourceExpr>
    >
{
    // cref if leaf, copy otherwise
    typename std::conditional<E1::is_leaf, const E1&, const E1>::type u_;

public:
    static constexpr bool is_leaf = false;

    //- Type to return for internal field
    typedef typename E1::value_type value_type;

    //- Type to return for containers
    typedef List_negate<typename E1::DiagExpr> DiagExpr;
    typedef List_negate<typename E1::UpperExpr> UpperExpr;
    typedef List_negate<typename E1::LowerExpr> LowerExpr;
    typedef List_negate<typename E1::FaceFluxExpr> FaceFluxExpr;
    typedef List_negate<typename E1::SourceExpr> SourceExpr;

    fvm_negate(const E1& u)
    :
        fvMatrixExpression
        <
            fvm_negate<E1>,
            DiagExpr,
            UpperExpr,
            LowerExpr,
            FaceFluxExpr,
            SourceExpr
        >
        (
            u.hasDiag(),
            u.hasUpper(),
            u.hasLower(),
            u.hasFaceFluxCorrection(),
            u.dimensions()
        ),
        u_(u)
    {}


    // lduMatrix part

        DiagExpr diag() const
        {
            return -u_.diag();
        }

        UpperExpr upper() const
        {
            return -u_.upper();
        }

        LowerExpr lower() const
        {
            return -u_.lower();
        }


    // fvMatrix

        SourceExpr source() const
        {
            return -u_.source();
        }

        UpperExpr internalCoeffs(const label i) const
        {
            return -u_.internalCoeffs(i);
        }

        UpperExpr boundaryCoeffs(const label i) const
        {
            return -u_.boundaryCoeffs(i);
        }

        FaceFluxExpr faceFluxCorrection() const
        {
            return -u_.faceFluxCorrection();
        }
};
template<typename E1>
fvm_negate<E1>
operator-
(
    fvMatrixExpression
    <
        E1,
        typename E1::DiagExpr,
        typename E1::UpperExpr,
        typename E1::LowerExpr,
        typename E1::FaceFluxExpr,
        typename E1::SourceExpr
    > const& u
)
{
    return fvm_negate<E1>(static_cast<const E1&>(u));
}


//- fvMatrix, internal parts and dimensions only. No upper/lower/boundaryCoeffs
//  etc. Used for source manipulation.
template<class E1, class E2>    //, class FaceFluxExpr>
class fvMatrixInternal
:
    public fvMatrixExpression
    <
        fvMatrixInternal<E1, E2>,
        E1, // Diag
        E1, // Upper, Not used
        E1, // Lower, Not used
        GeometricField<typename E1::value_type, fvPatchField, volMesh>, // FaceFlux, Not used
        E2  // Source
    >
{
public:

    static constexpr bool is_leaf = false;  //true;

    //- Type to return for internal field
    typedef typename E1::value_type value_type;

    //- Type to return for containers
    typedef E1 DiagExpr;
    typedef E1 UpperExpr;   // not yet used
    typedef E1 LowerExpr;   // not yet used
    typedef GeometricField<value_type, fvPatchField, volMesh>
        FaceFluxExpr;   // not yet used
    typedef E2 SourceExpr;


private:

    const E1 diag_;
    const SourceExpr source_;


public:

    // Construct from component expressions
    fvMatrixInternal
    (
        const dimensionSet& dimensions,
        const DiagExpr& diag,
        const SourceExpr& source
    )
    :
        fvMatrixExpression
        <
            fvMatrixInternal<DiagExpr, SourceExpr>,  //, FaceFluxExpr>,
            DiagExpr,
            UpperExpr,       // upper not yet used
            LowerExpr,       // lower not yet used
            FaceFluxExpr,    // face flux not used
            SourceExpr
        >
        (
            true,       //elems.hasDiag(),
            false,      //elems.hasUpper(),
            false,      //elems.hasLower(),
            false,      //elems.hasFaceFluxCorrection(),
            dimensions  //elems.dimensions()
        ),
        diag_(diag),
        source_(source)
    {
        //Pout<< "fvMatrixInternal :" << nl
        //    << "    diag_  : " << boost::core::demangle(typeid(diag_).name())
        //    << nl
        //    << "    source_: "
        //    << boost::core::demangle(typeid(source_).name())
        //    << nl
        //    << "    *this  : " << boost::core::demangle(typeid(*this).name())
        //    << endl;
    }

    // Construct from component expressions
    fvMatrixInternal
    (
        const dimensionSet& dimensions,
        const DiagExpr&& diag,
        const SourceExpr&& source
    )
    :
        fvMatrixExpression
        <
            fvMatrixInternal<DiagExpr, SourceExpr>,  //, FaceFluxExpr>,
            DiagExpr,
            UpperExpr,       // upper not yet used
            LowerExpr,       // lower not yet used
            FaceFluxExpr,       // faceflux not yet used
            SourceExpr
        >
        (
            true,       //elems.hasDiag(),
            false,      //elems.hasUpper(),
            false,      //elems.hasLower(),
            false,      //elems.hasFaceFluxCorrection(),
            dimensions  //elems.dimensions()
        ),
        diag_(diag),
        source_(source)
    {}


    // lduMatrix part

        auto diag() const
        {
            return diag_;
        }

        auto upper() const
        {
            return List<value_type>::null();
        }

        auto lower() const
        {
            return List<value_type>::null();
        }


    // fvMatrix

        auto source() const
        {
            return source_;
        }

        auto internalCoeffs(const label i) const
        {
            return List<value_type>::null();
        }

        auto boundaryCoeffs(const label i) const
        {
            return List<value_type>::null();
        }

        auto faceFluxCorrection() const
        {
            return GeometricField<value_type, fvPatchField, volMesh>::null();
        }
};
template<typename E1>
auto SuSp
(
    const dimensionSet& dimensions, // dimensions of expression
    const E1& expression,
    const GeometricField<typename E1::value_type, fvPatchField, volMesh>& fld
)
{
    //- Wrap of constant as a list expression
    typedef Expression::UniformListWrap<scalar> constant;
    //- Wrap of List as an expression
    typedef ListConstRefWrap<typename E1::value_type> expr;

    // Wrap zero
    const constant constantZero(fld.size(), 0.0);

    // Wrap mesh volume
    const expr V(fld.mesh().V());

    return fvMatrixInternal
    (
        dimVol*dimensions*fld.dimensions(),
        V*max(expression, constantZero),
        V*min(expression, constantZero)*expr(fld.internalField())
    );
}
template<typename E1>
auto Sp
(
    const dimensionSet& dimensions, // dimensions of expression
    const E1& expression,
    const GeometricField<typename E1::value_type, fvPatchField, volMesh>& fld
)
{
    //- Wrap of constant as a list expression
    typedef Expression::UniformListWrap<scalar> constant;
    //- Wrap of List as an expression
    typedef ListConstRefWrap<typename E1::value_type> expr;

    // Wrap zero
    const constant constantZero(fld.size(), 0.0);

    // Wrap mesh volume
    const expr V(fld.mesh().V());

    return fvMatrixInternal
    (
        dimVol*dimensions*fld.dimensions(),
        V*expr(fld.internalField()),
        constantZero    // Still provides an expression for source. TBD.
    );
}
template<typename E1>
auto Su
(
    const dimensionSet& dimensions, // dimensions of expression
    const E1& expression,
    const GeometricField<typename E1::value_type, fvPatchField, volMesh>& fld
)
{
    //- Wrap of constant as a list expression
    typedef Expression::UniformListWrap<scalar> constant;
    //- Wrap of List as an expression
    typedef ListConstRefWrap<typename E1::value_type> expr;

    // Wrap zero
    const constant constantZero(fld.size(), 0.0);

    // Wrap mesh volume
    const expr V(fld.mesh().V());

    return fvMatrixInternal
    (
        dimVol*dimensions*fld.dimensions(),
        constantZero    // Still provides an expression for diag. TBD.
        -V*expression
    );
}


//XXXXXXXX
// Alternative in-place fvMatrix sources. Can be removed once ones above
// work.
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
//XXXXXXXX



// Some discretisation using expressions
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template<typename E1, typename E2>
class interpolate
:
    public GeometricFieldExpression
    <
        interpolate<E1, E2>,

        //- Type to return for internal field
        List_interpolate
        <
            IndirectConstWrap<typename E1::IntExpr>,
            IndirectConstWrap<typename E1::IntExpr>,
            typename E2::IntExpr
        >,

        //- Type to return for uncoupled patch field
        typename E1::UncoupledPatchExpr,

        //- Type to return for coupled patch field
        List_interpolate
        <
            IndirectConstWrap<typename E1::IntExpr>,
            typename E1::CoupledPatchExpr,
            typename E2::CoupledPatchExpr
        >,

        //- Type of actual values
        typename E1::value_type
    >
{
    //- Internal values
    typename std::conditional
    <
        E1::is_leaf,
        const E1&,
        const E1
    >::type cellVals_;

    //- Interpolation weights
    typename std::conditional
    <
        E2::is_leaf,
        const E2&,
        const E2
    >::type faceWeights_;

    //- Addressing container
    const fvMesh& mesh_;

public:

    static constexpr bool is_leaf = false;

    //- Type to return for internal field
    typedef List_interpolate
    <
        IndirectConstWrap<typename E1::IntExpr>,
        IndirectConstWrap<typename E1::IntExpr>,
        typename E2::IntExpr
    > IntExpr;

    //- Type to return for uncoupled patch field
    typedef typename E1::UncoupledPatchExpr UncoupledPatchExpr;

    //- Type to return for coupled patch field
    typedef List_interpolate
    <
        IndirectConstWrap<typename E1::IntExpr>,
        typename E1::CoupledPatchExpr,
        typename E2::CoupledPatchExpr
    > CoupledPatchExpr;

    //- Type of actual values
    typedef typename E1::value_type value_type;


    // Construct from components
    interpolate
    (
        const E1& cellVals,
        const E2& faceWeights,
        const fvMesh& mesh
    )
    :
        GeometricFieldExpression
        <
            interpolate<E1, E2>,
            IntExpr,
            UncoupledPatchExpr,
            CoupledPatchExpr,
            value_type
        >
        (
            cellVals.dimensions()*faceWeights.dimensions(),
            faceWeights.oriented()
        ),
        cellVals_(cellVals),
        faceWeights_(faceWeights),
        mesh_(mesh)
    {}

    value_type operator[](const label i) const
    {
        const auto ownVal = cellVals_[mesh_.owner()[i]];
        const auto neiVal = cellVals_[mesh_.neighbour()[i]];
        return faceWeights_[i]*(ownVal-neiVal) + neiVal;
    }

    label size() const
    {
        return faceWeights_.size(); // number of internal faces
    }

    IntExpr internalField() const
    {
        return IntExpr
        (
            IndirectConstWrap<typename E1::IntExpr>
            (
                cellVals_.internalField(),
                mesh_.owner()
            ),
            IndirectConstWrap<typename E1::IntExpr>
            (
                cellVals_.internalField(),
                mesh_.neighbour()
            ),
            faceWeights_.internalField()
        );
    }

    UncoupledPatchExpr patchField(const label i) const
    {
        return UncoupledPatchExpr(cellVals_.patchField(i));
    }

    CoupledPatchExpr coupledPatchField(const label i) const
    {
        return CoupledPatchExpr
        (
            IndirectConstWrap<typename E1::IntExpr>
            (
                cellVals_.internalField(),
                mesh_.boundaryMesh()[i].faceCells()
            ),
            cellVals_.access
            (
                [&](const auto& fld, const label i)
                {
                    return ListConstTmpWrap<Field<value_type>>
                    (
                        fld.boundaryField()[i].patchNeighbourField()
                    );
                },
                i
            ),
            faceWeights_.coupledPatchField(i)
        );
    }
};


template<typename E1>
class linearInterpolate
:
    public interpolate
    <
        E1,
        GeometricFieldConstRefWrap<surfaceScalarField>
    >
{
public:

    // Construct from components
    linearInterpolate
    (
        const E1& cellVals,
        const fvMesh& mesh
    )
    :
        interpolate
        <
            E1,
            GeometricFieldConstRefWrap<surfaceScalarField>
        >
        (
            cellVals,
            mesh.surfaceInterpolation::weights().expr(),
            mesh
        )
    {}
};


template<class Type, class E1, class E2>
void fvmLaplacianUncorrected
(
    fvMatrix<Type>& fvm,
    const E1& gammaMagSf,
    const E2& deltaCoeffs
)
{
    typedef Expression::ListConstRefWrap<Type> expr;
    typedef Expression::ListRefWrap<Type> evaluator;

    const auto& vf = fvm.psi();

    // Set dimensions. Or should we just check?
    const_cast<dimensionSet&>(fvm.dimensions()).reset
    (
        deltaCoeffs.dimensions()
      * gammaMagSf.dimensions()
      * vf.dimensions()
    );

    Expression::ListRefWrap<scalar>
    (
        fvm.upper(),
        deltaCoeffs.internalField()*gammaMagSf.internalField()
    );
    fvm.negSumDiag();

    forAll(vf.boundaryField(), patchi)
    {
        const auto& pvf = vf.boundaryField()[patchi];
        auto& intCoeffs = fvm.internalCoeffs()[patchi];
        auto& bouCoeffs = fvm.boundaryCoeffs()[patchi];

        if (pvf.coupled())
        {
            // Evaluate deltaCoeffs
            scalarField pDeltaCoeffs;
            ListRefWrap<scalar>(pDeltaCoeffs, deltaCoeffs.patchField(patchi));

            //intCoeffs = pGamma*pvf.gradientInternalCoeffs(pDeltaCoeffs);
            evaluator
            (
                intCoeffs,
                expr(pvf.gradientInternalCoeffs(pDeltaCoeffs)())
               *gammaMagSf.patchField(patchi)
            );
            //bouCoeffs = -pGamma*pvf.gradientBoundaryCoeffs(pDeltaCoeffs);
            evaluator
            (
                bouCoeffs,
               -expr(pvf.gradientBoundaryCoeffs(pDeltaCoeffs)())
               *gammaMagSf.patchField(patchi)
            );
        }
        else
        {
            //intCoeffs = pGamma*pvf.gradientInternalCoeffs();
            evaluator
            (
                intCoeffs,
                expr(pvf.gradientInternalCoeffs()())
               *gammaMagSf.patchField(patchi)
            );
            //bouCoeffs = -pGamma*pvf.gradientBoundaryCoeffs();
            evaluator
            (
                bouCoeffs,
               -expr(pvf.gradientBoundaryCoeffs()())
               *gammaMagSf.patchField(patchi)
            );
        }
    }
}


} // End namespace Expression

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
