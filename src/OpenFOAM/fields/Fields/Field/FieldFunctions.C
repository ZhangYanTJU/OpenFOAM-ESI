/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "PstreamReduceOps.H"
#include "FieldReuseFunctions.H"

#define TEMPLATE template<class Type>
#include "FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * */

template<class Type>
void component
(
    Field<typename Field<Type>::cmptType>& result,
    const UList<Type>& f1,
    const direction d
)
{
    if (result.cdata_bytes() == f1.cdata_bytes())
    {
        // std::for_each
        TSEQ_FORALL_F_OP_F_FUNC_S_inplace(result, =, f1, .component, d)
    }
    else
    {
        // std::transform
        TSEQ_FORALL_F_OP_F_FUNC_S(result, =, f1, .component, d)
    }
}


template<class Type>
void T(Field<Type>& result, const UList<Type>& f1)
{
    if (result.cdata() == f1.cdata())
    {
        // std::for_each
        TSEQ_FORALL_F_OP_F_FUNC_inplace(result, =, f1, T)
    }
    else
    {
        // std::transform
        TSEQ_FORALL_F_OP_F_FUNC(result, =, f1, T)
    }
}


template<class Type, direction r>
void pow
(
    Field<typename powProduct<Type, r>::type>& result,
    const UList<Type>& f1
)
{
    typedef typename powProduct<Type, r>::type resultType;

    if (result.cdata_bytes() == f1.cdata_bytes())
    {
        // std::for_each
        TSEQ_FORALL_F_OP_FUNC_F_S_inplace
        (result, =, pow, f1, pTraits<resultType>::zero)
    }
    else
    {
        // std::transform
        TSEQ_FORALL_F_OP_FUNC_F_S
        (result, =, pow, f1, pTraits<resultType>::zero)
    }

}

template<class Type, direction r>
tmp<Field<typename powProduct<Type, r>::type>>
pow
(
    const UList<Type>& f1,
    typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type resultType;
    auto tres = tmp<Field<resultType>>::New(f1.size());
    pow<Type, r>(tres.ref(), f1);
    return tres;
}

template<class Type, direction r>
tmp<Field<typename powProduct<Type, r>::type>>
pow
(
    const tmp<Field<Type>>& tf1,
    typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type resultType;
    auto tres = reuseTmp<resultType, Type>::New(tf1);
    pow<Type, r>(tres.ref(), tf1());
    tf1.clear();
    return tres;
}


template<class Type>
void sqr
(
    Field<typename outerProduct<Type, Type>::type>& result,
    const UList<Type>& f1
)
{
    if (result.cdata_bytes() == f1.cdata_bytes())
    {
        // std::for_each
        TSEQ_FORALL_F_OP_FUNC_F_inplace(result, =, sqr, f1)
    }
    else
    {
        // std::transform
        TSEQ_FORALL_F_OP_FUNC_F(result, =, sqr, f1)
    }
}

template<class Type>
tmp<Field<typename outerProduct<Type, Type>::type>>
sqr(const UList<Type>& f1)
{
    typedef typename outerProduct<Type, Type>::type resultType;
    auto tres = tmp<Field<resultType>>::New(f1.size());
    sqr(tres.ref(), f1);
    return tres;
}

template<class Type>
tmp<Field<typename outerProduct<Type, Type>::type>>
sqr(const tmp<Field<Type>>& tf1)
{
    typedef typename outerProduct<Type, Type>::type resultType;
    auto tres = reuseTmp<resultType, Type>::New(tf1);
    sqr(tres.ref(), tf1());
    tf1.clear();
    return tres;
}


template<class Type>
void magSqr
(
    Field<typename typeOfMag<Type>::type>& result,
    const UList<Type>& f1
)
{
    if (result.cdata_bytes() == f1.cdata_bytes())
    {
        // std::for_each
        TSEQ_FORALL_F_OP_FUNC_F_inplace(result, =, magSqr, f1)
    }
    else
    {
        // std::transform
        TSEQ_FORALL_F_OP_FUNC_F(result, =, magSqr, f1)
    }
}

template<class Type>
tmp<Field<typename typeOfMag<Type>::type>>
magSqr(const UList<Type>& f1)
{
    typedef typename typeOfMag<Type>::type resultType;

    auto tres = tmp<Field<resultType>>::New(f1.size());
    magSqr(tres.ref(), f1);
    return tres;
}

template<class Type>
tmp<Field<typename typeOfMag<Type>::type>>
magSqr(const tmp<Field<Type>>& tf1)
{
    typedef typename typeOfMag<Type>::type resultType;

    auto tres = reuseTmp<resultType, Type>::New(tf1);
    magSqr(tres.ref(), tf1());
    tf1.clear();
    return tres;
}


template<class Type>
void mag
(
    Field<typename typeOfMag<Type>::type>& result,
    const UList<Type>& f1
)
{
    if (result.cdata_bytes() == f1.cdata_bytes())
    {
        // std::for_each
        TSEQ_FORALL_F_OP_FUNC_F_inplace(result, =, mag, f1)
    }
    else
    {
        // std::transform
        TSEQ_FORALL_F_OP_FUNC_F(result, =, mag, f1)
    }
}

template<class Type>
tmp<Field<typename typeOfMag<Type>::type>>
mag(const UList<Type>& f1)
{
    typedef typename typeOfMag<Type>::type resultType;

    auto tres = tmp<Field<resultType>>::New(f1.size());
    mag(tres.ref(), f1);
    return tres;
}

template<class Type>
tmp<Field<typename typeOfMag<Type>::type>>
mag(const tmp<Field<Type>>& tf1)
{
    typedef typename typeOfMag<Type>::type resultType;

    auto tres = reuseTmp<resultType, Type>::New(tf1);
    mag(tres.ref(), tf1());
    tf1.clear();
    return tres;
}


template<class Type>
void cmptMax
(
    Field<typename Field<Type>::cmptType>& result,
    const UList<Type>& f1
)
{
    if (result.cdata_bytes() == f1.cdata_bytes())
    {
        // std::for_each
        TSEQ_FORALL_F_OP_FUNC_F_inplace(result, =, cmptMax, f1)
    }
    else
    {
        // std::transform
        TSEQ_FORALL_F_OP_FUNC_F(result, =, cmptMax, f1)
    }
}

template<class Type>
tmp<Field<typename Field<Type>::cmptType>> cmptMax(const UList<Type>& f1)
{
    typedef typename Field<Type>::cmptType resultType;
    auto tres = tmp<Field<resultType>>::New(f1.size());
    cmptMax(tres.ref(), f1);
    return tres;
}

template<class Type>
tmp<Field<typename Field<Type>::cmptType>> cmptMax(const tmp<Field<Type>>& tf1)
{
    typedef typename Field<Type>::cmptType resultType;
    auto tres = reuseTmp<resultType, Type>::New(tf1);
    cmptMax(tres.ref(), tf1());
    tf1.clear();
    return tres;
}


template<class Type>
void cmptMin
(
    Field<typename Field<Type>::cmptType>& result,
    const UList<Type>& f1
)
{
    if (result.cdata_bytes() == f1.cdata_bytes())
    {
        // std::for_each
        TSEQ_FORALL_F_OP_FUNC_F_inplace(result, =, cmptMin, f1)
    }
    else
    {
        // std::transform
        TSEQ_FORALL_F_OP_FUNC_F(result, =, cmptMin, f1)
    }
}

template<class Type>
tmp<Field<typename Field<Type>::cmptType>> cmptMin(const UList<Type>& f1)
{
    typedef typename Field<Type>::cmptType resultType;
    auto tres = tmp<Field<resultType>>::New(f1.size());
    cmptMin(tres.ref(), f1);
    return tres;
}

template<class Type>
tmp<Field<typename Field<Type>::cmptType>> cmptMin(const tmp<Field<Type>>& tf1)
{
    typedef typename Field<Type>::cmptType resultType;
    auto tres = reuseTmp<resultType, Type>::New(tf1);
    cmptMin(tres.ref(), tf1());
    tf1.clear();
    return tres;
}


template<class Type>
void cmptAv
(
    Field<typename Field<Type>::cmptType>& result,
    const UList<Type>& f1
)
{
    if (result.cdata_bytes() == f1.cdata_bytes())
    {
        // std::for_each
        TSEQ_FORALL_F_OP_FUNC_F_inplace(result, =, cmptAv, f1)
    }
    else
    {
        // std::transform
        TSEQ_FORALL_F_OP_FUNC_F(result, =, cmptAv, f1)
    }
}

template<class Type>
tmp<Field<typename Field<Type>::cmptType>> cmptAv(const UList<Type>& f1)
{
    typedef typename Field<Type>::cmptType resultType;
    auto tres = tmp<Field<resultType>>::New(f1.size());
    cmptAv(tres.ref(), f1);
    return tres;
}

template<class Type>
tmp<Field<typename Field<Type>::cmptType>> cmptAv(const tmp<Field<Type>>& tf1)
{
    typedef typename Field<Type>::cmptType resultType;
    auto tres = reuseTmp<resultType, Type>::New(tf1);
    cmptAv(tres.ref(), tf1());
    tf1.clear();
    return tres;
}


template<class Type>
void cmptMag(Field<Type>& result, const UList<Type>& f1)
{
    if (result.cdata_bytes() == f1.cdata_bytes())
    {
        // std::for_each
        TSEQ_FORALL_F_OP_FUNC_F_inplace(result, =, cmptMag, f1)
    }
    else
    {
        // std::transform
        TSEQ_FORALL_F_OP_FUNC_F(result, =, cmptMag, f1)
    }
}

template<class Type>
tmp<Field<Type>> cmptMag(const UList<Type>& f1)
{
    auto tres = tmp<Field<Type>>::New(f1.size());
    cmptMag(tres.ref(), f1);
    return tres;
}

template<class Type>
tmp<Field<Type>> cmptMag(const tmp<Field<Type>>& tf1)
{
    auto tres = New(tf1);
    cmptMag(tres.ref(), tf1());
    tf1.clear();
    return tres;
}


template<class Type>
void cmptMagSqr(Field<Type>& result, const UList<Type>& f1)
{
    if (result.cdata_bytes() == f1.cdata_bytes())
    {
        // std::for_each
        TSEQ_FORALL_F_OP_FUNC_F_inplace(result, =, cmptMagSqr, f1)
    }
    else
    {
        // std::transform
        TSEQ_FORALL_F_OP_FUNC_F(result, =, cmptMagSqr, f1)
    }
}

template<class Type>
tmp<Field<Type>> cmptMagSqr(const UList<Type>& f1)
{
    auto tres = tmp<Field<Type>>::New(f1.size());
    cmptMagSqr(tres.ref(), f1);
    return tres;
}

template<class Type>
tmp<Field<Type>> cmptMagSqr(const tmp<Field<Type>>& tf1)
{
    auto tres = New(tf1);
    cmptMagSqr(tres.ref(), tf1());
    tf1.clear();
    return tres;
}


#define TMP_UNARY_FUNCTION(ReturnType, Func)                                   \
                                                                               \
template<class Type>                                                           \
ReturnType Func(const tmp<Field<Type>>& tf1)                                   \
{                                                                              \
    ReturnType res = Func(tf1());                                              \
    tf1.clear();                                                               \
    return res;                                                                \
}

template<class Type>
Type max(const UList<Type>& f1)
{
    if (f1.size())
    {
        Type result(f1[0]);
        TSEQ_FORALL_S_OP_FUNC_F_S(result, =, max, f1, result)
        return result;
    }

    return pTraits<Type>::min;
}

TMP_UNARY_FUNCTION(Type, max)

template<class Type>
Type min(const UList<Type>& f1)
{
    if (f1.size())
    {
        Type result(f1[0]);
        TSEQ_FORALL_S_OP_FUNC_F_S(result, =, min, f1, result)
        return result;
    }

    return pTraits<Type>::max;
}

TMP_UNARY_FUNCTION(Type, min)

template<class Type>
Type sum(const UList<Type>& f1)
{
    typedef typename Foam::typeOfSolve<Type>::type resultType;

    resultType result = Zero;

    if (f1.size())
    {
        // Use resultType() as functional cast
        TSEQ_FORALL_S_OP_FUNC_F(result, +=, resultType, f1)
    }

    return Type(result);
}

TMP_UNARY_FUNCTION(Type, sum)


// From MinMaxOps.H:
//   - Foam::minMax(const UList<Type>&)
//   - Foam::minMaxMag(const UList<Type>&)

TMP_UNARY_FUNCTION(MinMax<Type>, minMax)
TMP_UNARY_FUNCTION(scalarMinMax, minMaxMag)


template<class Type>
Type maxMagSqr(const UList<Type>& f1)
{
    if (f1.size())
    {
        Type result(f1[0]);
        TSEQ_FORALL_S_OP_FUNC_F_S
        (
            result,
            =,
            maxMagSqrOp<Type>(),
            f1,
            result
        )
        return result;
    }

    return Zero;
}

TMP_UNARY_FUNCTION(Type, maxMagSqr)

template<class Type>
Type minMagSqr(const UList<Type>& f1)
{
    if (f1.size())
    {
        Type result(f1[0]);
        TSEQ_FORALL_S_OP_FUNC_F_S
        (
            result,
            =,
            minMagSqrOp<Type>(),
            f1,
            result
        )
        return result;
    }

    return pTraits<Type>::rootMax;
}

TMP_UNARY_FUNCTION(Type, minMagSqr)

template<class Type>
typename scalarProduct<Type, Type>::type
sumProd(const UList<Type>& f1, const UList<Type>& f2)
{
    typedef typename scalarProduct<Type, Type>::type resultType;

    resultType result = Zero;
    if (f1.size() && (f1.size() == f2.size()))
    {
        // std::transform
        TSEQ_FORALL_S_OP_F_OP_F(result, +=, f1, &&, f2)
    }
    return result;
}


template<class Type>
Type sumCmptProd(const UList<Type>& f1, const UList<Type>& f2)
{
    Type result = Zero;
    if (f1.size() && (f1.size() == f2.size()))
    {
        TSEQ_FORALL_S_OP_FUNC_F_F
        (
            result,
            +=,
            cmptMultiply,
            f1,
            f2
        )
    }
    return result;
}


template<class Type>
typename outerProduct1<Type>::type
sumSqr(const UList<Type>& f1)
{
    typedef typename outerProduct1<Type>::type resultType;

    resultType result = Zero;
    if (f1.size())
    {
        TSEQ_FORALL_S_OP_FUNC_F(result, +=, sqr, f1)
    }
    return result;
}

template<class Type>
typename outerProduct1<Type>::type
sumSqr(const tmp<Field<Type>>& tf1)
{
    typedef typename outerProduct1<Type>::type resultType;
    resultType result = sumSqr(tf1());
    tf1.clear();
    return result;
}


template<class Type>
typename typeOfMag<Type>::type
sumMag(const UList<Type>& f1)
{
    typedef typename typeOfMag<Type>::type resultType;

    resultType result = Zero;
    if (f1.size())
    {
        TSEQ_FORALL_S_OP_FUNC_F(result, +=, mag, f1)
    }
    return result;
}

TMP_UNARY_FUNCTION(typename typeOfMag<Type>::type, sumMag)


template<class Type>
Type sumCmptMag(const UList<Type>& f1)
{
    Type result = Zero;
    if (f1.size())
    {
        TSEQ_FORALL_S_OP_FUNC_F(result, +=, cmptMag, f1)
    }
    return result;
}

TMP_UNARY_FUNCTION(Type, sumCmptMag)

template<class Type>
Type average(const UList<Type>& f1)
{
    if (f1.size())
    {
        Type result = sum(f1)/f1.size();

        return result;
    }

    WarningInFunction
        << "empty field, returning zero" << endl;

    return Zero;
}

TMP_UNARY_FUNCTION(Type, average)


// With reduction on ReturnType
#define G_UNARY_FUNCTION(ReturnType, gFunc, Func, rFunc)                       \
                                                                               \
template<class Type>                                                           \
ReturnType gFunc(const UList<Type>& f, const label comm)                       \
{                                                                              \
    ReturnType res = Func(f);                                                  \
    reduce(res, rFunc##Op<ReturnType>(), UPstream::msgType(), comm);           \
    return res;                                                                \
}                                                                              \
TMP_UNARY_FUNCTION(ReturnType, gFunc)

G_UNARY_FUNCTION(Type, gMax, max, max)
G_UNARY_FUNCTION(Type, gMin, min, min)
G_UNARY_FUNCTION(Type, gSum, sum, sum)
G_UNARY_FUNCTION(Type, gMaxMagSqr, maxMagSqr, maxMagSqr)
G_UNARY_FUNCTION(Type, gMinMagSqr, minMagSqr, minMagSqr)
G_UNARY_FUNCTION(Type, gSumCmptMag, sumCmptMag, sum)

G_UNARY_FUNCTION(MinMax<Type>, gMinMax, minMax, sum)
G_UNARY_FUNCTION(scalarMinMax, gMinMaxMag, minMaxMag, sum)

G_UNARY_FUNCTION(typename outerProduct1<Type>::type, gSumSqr, sumSqr, sum)
G_UNARY_FUNCTION(typename typeOfMag<Type>::type, gSumMag, sumMag, sum)

#undef G_UNARY_FUNCTION


template<class Type>
typename scalarProduct<Type, Type>::type gSumProd
(
    const UList<Type>& f1,
    const UList<Type>& f2,
    const label comm
)
{
    typedef typename scalarProduct<Type, Type>::type resultType;

    resultType result = sumProd(f1, f2);
    reduce(result, sumOp<resultType>(), UPstream::msgType(), comm);
    return result;
}

template<class Type>
Type gSumCmptProd
(
    const UList<Type>& f1,
    const UList<Type>& f2,
    const label comm
)
{
    Type result = sumCmptProd(f1, f2);
    reduce(result, sumOp<Type>(), UPstream::msgType(), comm);
    return result;
}

template<class Type>
Type gAverage
(
    const UList<Type>& f1,
    const label comm
)
{
    label n = f1.size();
    Type s = sum(f1);
    sumReduce(s, n, UPstream::msgType(), comm);

    if (n > 0)
    {
        Type result = s/n;

        return result;
    }

    WarningInFunction
        << "empty field, returning zero." << endl;

    return Zero;
}

TMP_UNARY_FUNCTION(Type, gAverage)

#undef TMP_UNARY_FUNCTION


// Implement BINARY_FUNCTION_TRANSFORM_FS for clamp
template<class Type>
void clamp
(
    Field<Type>& result,
    const UList<Type>& f1,
    const MinMax<Type>& range
)
{
    // Note: no checks for bad/invalid clamping ranges

    if (result.cdata() == f1.cdata())
    {
        // Apply in-place
        result.clamp_range(range);
    }
    else
    {
        std::transform
        (
            f1.cbegin(),
            f1.cbegin(result.size()),
            result.begin(),
            clampOp<Type>(range)
        );
    }
}

template<class Type>
void clamp
(
    Field<Type>& result,
    const UList<Type>& f1,
    const Foam::zero_one&   // Note: macros generate a const reference
)
{
    if (result.cdata() == f1.cdata())
    {
        // Apply in-place
        result.clamp_range(Foam::zero_one{});
    }
    else
    {
        std::transform
        (
            f1.cbegin(),
            f1.cbegin(result.size()),
            result.begin(),
            clampOp<Type>(Foam::zero_one{})
        );
    }
}

BINARY_FUNCTION_INTERFACE_FS(Type, Type, MinMax<Type>, clamp)
BINARY_FUNCTION_INTERFACE_FS(Type, Type, Foam::zero_one, clamp)


BINARY_FUNCTION(Type, Type, Type, max)
BINARY_FUNCTION(Type, Type, Type, min)
BINARY_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_FUNCTION(Type, Type, Type, cmptDivide)

BINARY_TYPE_FUNCTION(Type, Type, Type, max)
BINARY_TYPE_FUNCTION(Type, Type, Type, min)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptDivide)

BINARY_TYPE_FUNCTION_FS(Type, Type, MinMax<Type>, clip)  // Same as clamp

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

TERNARY_FUNCTION(Type, Type, Type, scalar, lerp)
TERNARY_TYPE_FUNCTION_FFS(Type, Type, Type, scalar, lerp)


/* * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * */

UNARY_OPERATOR(Type, Type, -, negate)

BINARY_OPERATOR(Type, Type, scalar, *, multiply)
BINARY_OPERATOR(Type, scalar, Type, *, multiply)
BINARY_OPERATOR(Type, Type, scalar, /, divide)

BINARY_TYPE_OPERATOR_SF(Type, scalar, Type, *, multiply)
BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, *, multiply)

BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, /, divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define PRODUCT_OPERATOR(product, Op, OpFunc)                                  \
                                                                               \
template<class Type1, class Type2>                                             \
void OpFunc                                                                    \
(                                                                              \
    Field<typename product<Type1, Type2>::type>& result,                       \
    const UList<Type1>& f1,                                                    \
    const UList<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    if                                                                         \
    (                                                                          \
        result.cdata_bytes() == f1.cdata_bytes()                               \
     || result.cdata_bytes() == f2.cdata_bytes()                               \
    )                                                                          \
    {                                                                          \
        /* std::for_each */                                                    \
        TSEQ_FORALL_F_OP_F_OP_F_inplace(result, =, f1, Op, f2)                 \
    }                                                                          \
    else                                                                       \
    {                                                                          \
        /* std::transform */                                                   \
        TSEQ_FORALL_F_OP_F_OP_F(result, =, f1, Op, f2)                         \
    }                                                                          \
}                                                                              \
                                                                               \
template<class Type1, class Type2>                                             \
tmp<Field<typename product<Type1, Type2>::type>>                               \
operator Op(const UList<Type1>& f1, const UList<Type2>& f2)                    \
{                                                                              \
    typedef typename product<Type1, Type2>::type resultType;                   \
    auto tres = tmp<Field<resultType>>::New(f1.size());                        \
    OpFunc(tres.ref(), f1, f2);                                                \
    return tres;                                                               \
}                                                                              \
                                                                               \
template<class Type1, class Type2>                                             \
tmp<Field<typename product<Type1, Type2>::type>>                               \
operator Op(const UList<Type1>& f1, const tmp<Field<Type2>>& tf2)              \
{                                                                              \
    typedef typename product<Type1, Type2>::type resultType;                   \
    auto tres = reuseTmp<resultType, Type2>::New(tf2);                         \
    OpFunc(tres.ref(), f1, tf2());                                             \
    tf2.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
template<class Type1, class Type2>                                             \
tmp<Field<typename product<Type1, Type2>::type>>                               \
operator Op(const tmp<Field<Type1>>& tf1, const UList<Type2>& f2)              \
{                                                                              \
    typedef typename product<Type1, Type2>::type resultType;                   \
    auto tres = reuseTmp<resultType, Type1>::New(tf1);                         \
    OpFunc(tres.ref(), tf1(), f2);                                             \
    tf1.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
template<class Type1, class Type2>                                             \
tmp<Field<typename product<Type1, Type2>::type>>                               \
operator Op(const tmp<Field<Type1>>& tf1, const tmp<Field<Type2>>& tf2)        \
{                                                                              \
    typedef typename product<Type1, Type2>::type resultType;                   \
    auto tres = reuseTmpTmp<resultType, Type1, Type1, Type2>::New(tf1, tf2);   \
    OpFunc(tres.ref(), tf1(), tf2());                                          \
    tf1.clear();                                                               \
    tf2.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
template<class Type, class Form, class Cmpt, direction nCmpt>                  \
void OpFunc                                                                    \
(                                                                              \
    Field<typename product<Type, Form>::type>& result,                         \
    const UList<Type>& f1,                                                     \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
)                                                                              \
{                                                                              \
    if                                                                         \
    (                                                                          \
        result.cdata_bytes() == f1.cdata_bytes()                               \
    )                                                                          \
    {                                                                          \
        /* std::for_each */                                                    \
        TSEQ_FORALL_F_OP_F_OP_S_inplace                                        \
        (result, =, f1, Op, static_cast<const Form&>(vs))                      \
    }                                                                          \
    else                                                                       \
    {                                                                          \
        /* std::transform */                                                   \
        TSEQ_FORALL_F_OP_F_OP_S                                                \
        (result, =, f1, Op, static_cast<const Form&>(vs))                      \
    }                                                                          \
}                                                                              \
                                                                               \
template<class Type, class Form, class Cmpt, direction nCmpt>                  \
tmp<Field<typename product<Type, Form>::type>>                                 \
operator Op(const UList<Type>& f1, const VectorSpace<Form,Cmpt,nCmpt>& vs)     \
{                                                                              \
    typedef typename product<Type, Form>::type resultType;                     \
    auto tres = tmp<Field<resultType>>::New(f1.size());                        \
    OpFunc(tres.ref(), f1, static_cast<const Form&>(vs));                      \
    return tres;                                                               \
}                                                                              \
                                                                               \
template<class Type, class Form, class Cmpt, direction nCmpt>                  \
tmp<Field<typename product<Type, Form>::type>>                                 \
operator Op                                                                    \
(                                                                              \
    const tmp<Field<Type>>& tf1,                                               \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type resultType;                     \
    auto tres = reuseTmp<resultType, Type>::New(tf1);                          \
    OpFunc(tres.ref(), tf1(), static_cast<const Form&>(vs));                   \
    tf1.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
template<class Form, class Cmpt, direction nCmpt, class Type>                  \
void OpFunc                                                                    \
(                                                                              \
    Field<typename product<Form, Type>::type>& result,                         \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const UList<Type>& f1                                                      \
)                                                                              \
{                                                                              \
    if                                                                         \
    (                                                                          \
        result.cdata_bytes() == f1.cdata_bytes()                               \
    )                                                                          \
    {                                                                          \
        /* std::for_each */                                                    \
        TSEQ_FORALL_F_OP_S_OP_F_inplace                                        \
        (result, =, static_cast<const Form&>(vs), Op, f1)                      \
    }                                                                          \
    else                                                                       \
    {                                                                          \
        /* std::transform */                                                   \
        TSEQ_FORALL_F_OP_S_OP_F                                                \
        (result, =, static_cast<const Form&>(vs), Op, f1)                      \
    }                                                                          \
}                                                                              \
                                                                               \
template<class Form, class Cmpt, direction nCmpt, class Type>                  \
tmp<Field<typename product<Form, Type>::type>>                                 \
operator Op(const VectorSpace<Form,Cmpt,nCmpt>& vs, const UList<Type>& f1)     \
{                                                                              \
    typedef typename product<Form, Type>::type resultType;                     \
    auto tres = tmp<Field<resultType>>::New(f1.size());                        \
    OpFunc(tres.ref(), static_cast<const Form&>(vs), f1);                      \
    return tres;                                                               \
}                                                                              \
                                                                               \
template<class Form, class Cmpt, direction nCmpt, class Type>                  \
tmp<Field<typename product<Form, Type>::type>>                                 \
operator Op                                                                    \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& vs, const tmp<Field<Type>>& tf1        \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type resultType;                     \
    auto tres = reuseTmp<resultType, Type>::New(tf1);                          \
    OpFunc(tres.ref(), static_cast<const Form&>(vs), tf1());                   \
    tf1.clear();                                                               \
    return tres;                                                               \
}

PRODUCT_OPERATOR(typeOfSum, +, add)
PRODUCT_OPERATOR(typeOfSum, -, subtract)

PRODUCT_OPERATOR(outerProduct, *, outer)
PRODUCT_OPERATOR(crossProduct, ^, cross)
PRODUCT_OPERATOR(innerProduct, &, dot)
PRODUCT_OPERATOR(scalarProduct, &&, dotdot)

#undef PRODUCT_OPERATOR


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
