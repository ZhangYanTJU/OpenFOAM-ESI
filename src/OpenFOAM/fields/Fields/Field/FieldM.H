/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2022-2024 OpenCFD Ltd.
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
    Declaration macros for Field\<Type\> algebra.

\*---------------------------------------------------------------------------*/

#ifndef Foam_FieldM_H
#define Foam_FieldM_H

#include "errorCheckFields.H"   // Field size checks (fulldebug mode)
#include "ListLoopM.H"          // List access macros

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ------------
// Unary Free Function : f1 OP Func(f2)

#define TSEQ_FORALL_F_OP_FUNC_F_impl(f1, OP, FUNC, f2, Restrict)               \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP " " #FUNC "(f2)");                           \
                                                                               \
    /* Field access */                                                         \
    TFOR_List_ACCESS(f1, f1P, Restrict);                                       \
    TFOR_List_CONST_ACCESS(f2, f2P, Restrict);                                 \
    TFOR_List_LENGTH(f1, loop_len);                                            \
                                                                               \
    /* Loop: f1 OP FUNC(f2) */                                                 \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP FUNC(f2P[i]);                                              \
    }                                                                          \
}

#define TSEQ_FORALL_F_OP_FUNC_F(f1, OP, FUNC, f2) \
        TSEQ_FORALL_F_OP_FUNC_F_impl(f1, OP, FUNC, f2, __restrict__)

#define TSEQ_FORALL_F_OP_FUNC_F_inplace(f1, OP, FUNC, f2) \
        TSEQ_FORALL_F_OP_FUNC_F_impl(f1, OP, FUNC, f2,)

// old-style
#define TFOR_ALL_F_OP_FUNC_F(typeF1, f1, OP, FUNC, typeF2, f2) \
     TSEQ_FORALL_F_OP_FUNC_F_impl(f1, OP, FUNC, f2, __restrict__)


// ------------
// Nullary Member Function : f1 OP f2.FUNC()
// NB: the '.' is not passed in the func parameter

#define TSEQ_FORALL_F_OP_F_FUNC_impl(f1, OP, f2, FUNC, Restrict)               \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP " f2 ." #FUNC "()");                         \
                                                                               \
    /* Field access */                                                         \
    TFOR_List_ACCESS(f1, f1P, Restrict);                                       \
    TFOR_List_CONST_ACCESS(f2, f2P, Restrict);                                 \
    TFOR_List_LENGTH(f1, loop_len);                                            \
                                                                               \
    /* Loop: f1 OP f2.FUNC() */                                                \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP (f2P[i]).FUNC();                                           \
    }                                                                          \
}

#define TSEQ_FORALL_F_OP_F_FUNC(f1, OP, f2, FUNC) \
        TSEQ_FORALL_F_OP_F_FUNC_impl(f1, OP, f2, FUNC, __restrict__)

#define TSEQ_FORALL_F_OP_F_FUNC_inplace(f1, OP, f2, FUNC) \
        TSEQ_FORALL_F_OP_F_FUNC_impl(f1, OP, f2, FUNC,)

// old-style
#define TFOR_ALL_F_OP_F_FUNC(typeF1, f1, OP, typeF2, f2, FUNC) \
     TSEQ_FORALL_F_OP_F_FUNC_impl(f1, OP, f2, FUNC, __restrict__)


// ------------
// Binary Free Function : f1 OP FUNC(f2, f3)

#define TSEQ_FORALL_F_OP_FUNC_F_F_impl(f1, OP, FUNC, f2, f3, Restrict)         \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, f3, "f1 " #OP " " #FUNC "(f2, f3)");                   \
                                                                               \
    /* Field access */                                                         \
    TFOR_List_ACCESS(f1, f1P, Restrict);                                       \
    TFOR_List_CONST_ACCESS(f2, f2P, Restrict);                                 \
    TFOR_List_CONST_ACCESS(f3, f3P, Restrict);                                 \
    TFOR_List_LENGTH(f1, loop_len);                                            \
                                                                               \
    /* Loop: f1 OP FUNC(f2, f3) */                                             \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP FUNC((f2P[i]), (f3P[i]));                                  \
    }                                                                          \
}

#define TSEQ_FORALL_F_OP_FUNC_F_F(f1, OP, FUNC, f2, f3) \
        TSEQ_FORALL_F_OP_FUNC_F_F_impl(f1, OP, FUNC, f2, f3, __restrict__)

#define TSEQ_FORALL_F_OP_FUNC_F_F_inplace(f1, OP, FUNC, f2, f3) \
        TSEQ_FORALL_F_OP_FUNC_F_F_impl(f1, OP, FUNC, f2, f3,)

// old-style
#define TFOR_ALL_F_OP_FUNC_F_F(typeF1, f1, OP, FUNC, typeF2, f2, typeF3, f3) \
     TSEQ_FORALL_F_OP_FUNC_F_F_impl(f1, OP, FUNC, f2, f3, __restrict__)


// ------------
// [reduction] Binary Free Function : s OP FUNC(f1, f2)

#define TSEQ_FORALL_S_OP_FUNC_F_F_impl(s, OP, FUNC, f1, f2, Unused)            \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "s " #OP " " #FUNC "(f1, f2)");                        \
                                                                               \
    /* Field access */                                                         \
    TFOR_List_CONST_ACCESS(f1, f1P, __restrict__);                             \
    TFOR_List_CONST_ACCESS(f2, f2P, __restrict__);                             \
    TFOR_List_LENGTH(f1, loop_len);                                            \
                                                                               \
    /* Loop: s OP FUNC(f1, f2) */                                              \
                                                                               \
    /* pragmas, reduction... */                                                \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (s) OP FUNC((f1P[i]), (f2P[i]));                                       \
    }                                                                          \
}

#define TSEQ_FORALL_S_OP_FUNC_F_F(s, OP, FUNC, f1, f2) \
        TSEQ_FORALL_S_OP_FUNC_F_F_impl(s, OP, FUNC, f1, f2, __restrict__)

#define TSEQ_FORALL_S_OP_FUNC_F_F_inplace(s, OP, FUNC, f1, f2) \
        TSEQ_FORALL_S_OP_FUNC_F_F_impl(s, OP, FUNC, f1, f2,)

// old-style
#define TFOR_ALL_S_OP_FUNC_F_F(typeS, s, OP, FUNC, typeF1, f1, typeF2, f2) \
     TSEQ_FORALL_S_OP_FUNC_F_F_impl(s, OP, FUNC, f1, f2, __restrict__)


// ------------
// Binary Free Function : f1 OP FUNC(f2, s)

#define TSEQ_FORALL_F_OP_FUNC_F_S_impl(f1, OP, FUNC, f2, s, Restrict)          \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP " " #FUNC "(f2, s)");                        \
                                                                               \
    /* Field access */                                                         \
    TFOR_List_ACCESS(f1, f1P, Restrict);                                       \
    TFOR_List_CONST_ACCESS(f2, f2P, Restrict);                                 \
    TFOR_List_LENGTH(f1, loop_len);                                            \
                                                                               \
    /* Loop: f1 OP FUNC(f2, s) */                                              \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP FUNC((f2P[i]), (s));                                       \
    }                                                                          \
}

#define TSEQ_FORALL_F_OP_FUNC_F_S(f1, OP, FUNC, f2, s) \
        TSEQ_FORALL_F_OP_FUNC_F_S_impl(f1, OP, FUNC, f2, s, __restrict__)

#define TSEQ_FORALL_F_OP_FUNC_F_S_inplace(f1, OP, FUNC, f2, s) \
        TSEQ_FORALL_F_OP_FUNC_F_S_impl(f1, OP, FUNC, f2, s,)

// old-style
#define TFOR_ALL_F_OP_FUNC_F_S(typeF1, f1, OP, FUNC, typeF2, f2, typeS, s) \
     TSEQ_FORALL_F_OP_FUNC_F_S_impl(f1, OP, FUNC, f2, s, __restrict__)


// ------------
// [reduction] Binary Free Function : s1 OP FUNC(f, s2)

#define TSEQ_FORALL_S_OP_FUNC_F_S_impl(s1, OP, FUNC, f, s2, Unused)            \
{                                                                              \
    /* Field access */                                                         \
    TFOR_List_CONST_ACCESS(f, fP, __restrict__);                               \
    TFOR_List_LENGTH(f, loop_len);                                             \
                                                                               \
    /* Loop: s1 OP FUNC(f, s2) */                                              \
                                                                               \
    /* pragmas, reduction... */                                                \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (s1) OP FUNC((fP[i]), (s2));                                           \
    }                                                                          \
}

#define TSEQ_FORALL_S_OP_FUNC_F_S(s1, OP, FUNC, f, s2) \
        TSEQ_FORALL_S_OP_FUNC_F_S_impl(s1, OP, FUNC, f, s2, __restrict__)

#define TSEQ_FORALL_S_OP_FUNC_F_S_inplace(s1, OP, FUNC, f, s2) \
        TSEQ_FORALL_S_OP_FUNC_F_S_impl(s1, OP, FUNC, f, s2,)

// old-style
#define TFOR_ALL_S_OP_FUNC_F_S(typeS1, s1, OP, FUNC, typeF, f, typeS2, s2) \
     TSEQ_FORALL_S_OP_FUNC_F_S_impl(s1, OP, FUNC, f, s2, __restrict__)


// ------------
// Binary Free Function : f1 OP FUNC(s, f2)

#define TSEQ_FORALL_F_OP_FUNC_S_F_impl(f1, OP, FUNC, s, f2, Restrict)          \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP " " #FUNC "(s, f2)");                        \
                                                                               \
    /* Field access */                                                         \
    TFOR_List_ACCESS(f1, f1P, Restrict);                                       \
    TFOR_List_CONST_ACCESS(f2, f2P, Restrict);                                 \
    TFOR_List_LENGTH(f1, loop_len);                                            \
                                                                               \
    /* Loop: f1 OP1 f2 OP2 f3 */                                               \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP FUNC((s), (f2P[i]));                                       \
    }                                                                          \
}

#define TSEQ_FORALL_F_OP_FUNC_S_F(f1, OP, FUNC, s, f2) \
        TSEQ_FORALL_F_OP_FUNC_S_F_impl(f1, OP, FUNC, s, f2, __restrict__)

#define TSEQ_FORALL_F_OP_FUNC_S_F_inplace(f1, OP, FUNC, s, f2) \
        TSEQ_FORALL_F_OP_FUNC_S_F_impl(f1, OP, FUNC, s, f2,)

// old-style
#define TFOR_ALL_F_OP_FUNC_S_F(typeF1, f1, OP, FUNC, typeS, s, typeF2, f2) \
     TSEQ_FORALL_F_OP_FUNC_S_F_impl(f1, OP, FUNC, s, f2, __restrict__)


// ------------
// Binary Free Function : f1 OP FUNC(s1, s2)

#define TSEQ_FORALL_F_OP_FUNC_S_S_impl(f1, OP, FUNC, s1, s2, Unused)           \
{                                                                              \
    /* Field access */                                                         \
    TFOR_List_ACCESS(f1, f1P, __restrict__);                                   \
    TFOR_List_LENGTH(f1, loop_len);                                            \
                                                                               \
    /* Loop: f1 OP FUNC(s1, s2) */                                             \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP FUNC((s1), (s2));                                          \
    }                                                                          \
}

#define TSEQ_FORALL_F_OP_FUNC_S_S(f1, OP, FUNC, s1, s2) \
        TSEQ_FORALL_F_OP_FUNC_S_S_impl(f1, OP, FUNC, s1, s2, __restrict__)

#define TSEQ_FORALL_F_OP_FUNC_S_S_inplace(f1, OP, FUNC, s1, s2) \
        TSEQ_FORALL_F_OP_FUNC_S_S_impl(f1, OP, FUNC, s1, s2,)

// old-style
#define TFOR_ALL_F_OP_FUNC_S_S(typeF1, f1, OP, FUNC, typeS1, s1, typeS2, s2) \
     TSEQ_FORALL_F_OP_FUNC_S_S_impl(f1, OP, FUNC, s1, s2, __restrict__)


// ------------
// Unary Member Function : f1 OP f2.FUNC(s)
// NB: the '.' is passed in the func parameter

#define TSEQ_FORALL_F_OP_F_FUNC_S_impl(f1, OP, f2, FUNC, s, Restrict)          \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP " f2 " #FUNC "(s)");                         \
                                                                               \
    /* Field access */                                                         \
    TFOR_List_ACCESS(f1, f1P, Restrict);                                       \
    TFOR_List_CONST_ACCESS(f2, f2P, Restrict);                                 \
    TFOR_List_LENGTH(f1, loop_len);                                            \
                                                                               \
    /* Loop: f1 OP f2 FUNC(s) */                                               \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP (f2P[i]) FUNC((s));                                        \
    }                                                                          \
}

#define TSEQ_FORALL_F_OP_F_FUNC_S(f1, OP, f2, FUNC, s) \
        TSEQ_FORALL_F_OP_F_FUNC_S_impl(f1, OP, f2, FUNC, s, __restrict__)

#define TSEQ_FORALL_F_OP_F_FUNC_S_inplace(f1, OP, f2, FUNC, s) \
        TSEQ_FORALL_F_OP_F_FUNC_S_impl(f1, OP, f2, FUNC, s,)

// old-style
#define TFOR_ALL_F_OP_F_FUNC_S(typeF1, f1, OP, typeF2, f2, FUNC, typeS, s) \
     TSEQ_FORALL_F_OP_F_FUNC_S_impl(f1, OP, f2, FUNC, s, __restrict__)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ------------
// Ternary Free Function : f1 OP FUNC(f2, f3, f4)

#define TSEQ_FORALL_F_OP_FUNC_F_F_F_impl(f1, OP, FUNC, f2, f3, f4, Restrict)   \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, f3, f4, "f1 " #OP " " #FUNC "(f2, f3, f4)");           \
                                                                               \
    /* Field access */                                                         \
    TFOR_List_ACCESS(f1, f1P, Restrict);                                       \
    TFOR_List_CONST_ACCESS(f2, f2P, Restrict);                                 \
    TFOR_List_CONST_ACCESS(f3, f3P, Restrict);                                 \
    TFOR_List_CONST_ACCESS(f4, f4P, Restrict);                                 \
    TFOR_List_LENGTH(f1, loop_len);                                            \
                                                                               \
    /* Loop: f1 OP FUNC(f2, f3, f4) */                                         \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP FUNC((f2P[i]), (f3P[i]), (f4P[i]));                        \
    }                                                                          \
}

#define TSEQ_FORALL_F_OP_FUNC_F_F_F(f1, OP, FUNC, f2, f3, f4) \
        TSEQ_FORALL_F_OP_FUNC_F_F_F_impl(f1, OP, FUNC, f2, f3, f4, __restrict__)

#define TSEQ_FORALL_F_OP_FUNC_F_F_F_inplace(f1, OP, FUNC, f2, f3, f4) \
        TSEQ_FORALL_F_OP_FUNC_F_F_F_impl(f1, OP, FUNC, f2, f3, f4,)

// old-style
#define TFOR_ALL_F_OP_FUNC_F_F_F(type1, f1, OP, FUNC, type2, f2, type3, f3, type4, f4) \
     TSEQ_FORALL_F_OP_FUNC_F_F_F_impl(f1, OP, FUNC, f2, f3, f4, __restrict__)


// ------------
// Ternary Free Function : f1 OP FUNC(f2, f3, s4)

#define TSEQ_FORALL_F_OP_FUNC_F_F_S_impl(f1, OP, FUNC, f2, f3, s4, Restrict)   \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, f3, "f1 " #OP " " #FUNC "(f2, f3, s)");                \
                                                                               \
    /* Field access */                                                         \
    TFOR_List_ACCESS(f1, f1P, Restrict);                                       \
    TFOR_List_CONST_ACCESS(f2, f2P, Restrict);                                 \
    TFOR_List_CONST_ACCESS(f3, f3P, Restrict);                                 \
    TFOR_List_LENGTH(f1, loop_len);                                            \
                                                                               \
    /* Loop: f1 OP FUNC(f2, f3, s4) */                                         \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP FUNC((f2P[i]), (f3P[i]), (s4));                            \
    }                                                                          \
}

#define TSEQ_FORALL_F_OP_FUNC_F_F_S(f1, OP, FUNC, f2, f3, s4) \
        TSEQ_FORALL_F_OP_FUNC_F_F_S_impl(f1, OP, FUNC, f2, f3, s4, __restrict__)

#define TSEQ_FORALL_F_OP_FUNC_F_F_S_inplace(f1, OP, FUNC, f2, f3, s4) \
        TSEQ_FORALL_F_OP_FUNC_F_F_S_impl(f1, OP, FUNC, f2, f3, s4,)

// old-style
#define TFOR_ALL_F_OP_FUNC_F_F_S(type1, f1, OP, FUNC, type2, f2, type3, f3, type4, s4) \
     TSEQ_FORALL_F_OP_FUNC_F_F_S_impl(f1, OP, FUNC, f2, f3, s4, __restrict__)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ------------
// Member operator : this field f1 OP1 f2 OP2 f3

#define TSEQ_FORALL_F_OP_F_OP_F_impl(f1, OP1, f2, OP2, f3, Restrict)           \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, f3, "f1 " #OP1 " f2 " #OP2 " f3");                     \
                                                                               \
    /* Field access */                                                         \
    TFOR_List_ACCESS(f1, f1P, Restrict);                                       \
    TFOR_List_CONST_ACCESS(f2, f2P, Restrict);                                 \
    TFOR_List_CONST_ACCESS(f3, f3P, Restrict);                                 \
    TFOR_List_LENGTH(f1, loop_len);                                            \
                                                                               \
    /* Loop: f1 OP1 f2 OP2 f3 */                                               \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP1 (f2P[i]) OP2 (f3P[i]);                                    \
    }                                                                          \
}

#define TSEQ_FORALL_F_OP_F_OP_F(f1, OP1, f2, OP2, f3) \
        TSEQ_FORALL_F_OP_F_OP_F_impl(f1, OP1, f2, OP2, f3, __restrict__)

#define TSEQ_FORALL_F_OP_F_OP_F_inplace(f1, OP1, f2, OP2, f3) \
        TSEQ_FORALL_F_OP_F_OP_F_impl(f1, OP1, f2, OP2, f3,)

// old-style
#define TFOR_ALL_F_OP_F_OP_F(type1, f1, OP1, type2, f2, OP2, type3, f3) \
     TSEQ_FORALL_F_OP_F_OP_F_impl(f1, OP1, f2, OP2, f3, __restrict__)


// ------------
// Member operator : this field f1 OP1 s OP2 f2

#define TSEQ_FORALL_F_OP_S_OP_F_impl(f1, OP1, s, OP2, f2, Restrict)            \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP1 " s " #OP2 " f2");                          \
                                                                               \
    /* Field access */                                                         \
    TFOR_List_ACCESS(f1, f1P, Restrict);                                       \
    TFOR_List_CONST_ACCESS(f2, f2P, Restrict);                                 \
    TFOR_List_LENGTH(f1, loop_len);                                            \
                                                                               \
    /* Loop: f1 OP1 s OP2 f2 */                                                \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP1 (s) OP2 (f2P[i]);                                         \
    }                                                                          \
}

#define TSEQ_FORALL_F_OP_S_OP_F(f1, OP1, s, OP2, f2) \
        TSEQ_FORALL_F_OP_S_OP_F_impl(f1, OP1, s, OP2, f2, __restrict__)

#define TSEQ_FORALL_F_OP_S_OP_F_inplace(f1, OP1, s, OP2, f2) \
        TSEQ_FORALL_F_OP_S_OP_F_impl(f1, OP1, s, OP2, f2,)

// old-style
#define TFOR_ALL_F_OP_S_OP_F(type1, f1, OP1, typeS, s, OP2, type2, f2) \
    TSEQ_FORALL_F_OP_S_OP_F_impl(f1, OP1, s, OP2, f2, __restrict__)


// ------------
// Member operator : this field f1 OP1 f2 OP2 s

#define TSEQ_FORALL_F_OP_F_OP_S_impl(f1, OP1, f2, OP2, s, Restrict)            \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP1 " f2 " #OP2 " s");                          \
                                                                               \
    /* Field access */                                                         \
    TFOR_List_ACCESS(f1, f1P, Restrict);                                       \
    TFOR_List_CONST_ACCESS(f2, f2P, Restrict);                                 \
    TFOR_List_LENGTH(f1, loop_len);                                            \
                                                                               \
    /* Loop f1 OP1 s OP2 f2 */                                                 \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP1 (f2P[i]) OP2 (s);                                         \
    }                                                                          \
}

#define TSEQ_FORALL_F_OP_F_OP_S(f1, OP1, f2, OP2, s) \
        TSEQ_FORALL_F_OP_F_OP_S_impl(f1, OP1, f2, OP2, s, __restrict__)

#define TSEQ_FORALL_F_OP_F_OP_S_inplace(f1, OP1, f2, OP2, s) \
        TSEQ_FORALL_F_OP_F_OP_S_impl(f1, OP1, f2, OP2, s,)

// old-style
#define TFOR_ALL_F_OP_F_OP_S(type1, f1, OP1, type2, f2, OP2, typeS, s) \
     TSEQ_FORALL_F_OP_F_OP_S_impl(f1, OP1, f2, OP2, s, __restrict__)


// ------------
// Member operator : this field f1 OP f2

#define TSEQ_FORALL_F_OP_F_impl(f1, OP, f2, Restrict)                          \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP " f2");                                      \
                                                                               \
    /* Field access */                                                         \
    TFOR_List_ACCESS(f1, f1P, Restrict);                                       \
    TFOR_List_CONST_ACCESS(f2, f2P, Restrict);                                 \
    TFOR_List_LENGTH(f1, loop_len);                                            \
                                                                               \
    /* Loop: f1 OP f2 */                                                       \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP (f2P[i]);                                                  \
    }                                                                          \
}

#define TSEQ_FORALL_F_OP_F(f1, OP, f2) \
        TSEQ_FORALL_F_OP_F_impl(f1, OP, f2, __restrict__)

#define TSEQ_FORALL_F_OP_F_inplace(f1, OP, f2) \
        TSEQ_FORALL_F_OP_F_impl(f1, OP, f2,)

// old-style
#define TFOR_ALL_F_OP_F(type1, f1, OP, type2, f2) \
     TSEQ_FORALL_F_OP_F_impl(f1, OP, f2, __restrict__)


// ------------
// Member operator : this field f1 OP1 OP2 f2

#define TSEQ_FORALL_F_OP_OP_F_impl(f1, OP1, OP2, f2, Restrict)                 \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, #OP1 " " #OP2 " f2");                                  \
                                                                               \
    /* Field access */                                                         \
    TFOR_List_ACCESS(f1, f1P, Restrict);                                       \
    TFOR_List_CONST_ACCESS(f2, f2P, Restrict);                                 \
    TFOR_List_LENGTH(f1, loop_len);                                            \
                                                                               \
    /* Loop: f1 OP1 OP2 f2 */                                                  \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP1 OP2 (f2P[i]);                                             \
    }                                                                          \
}

#define TSEQ_FORALL_F_OP_OP_F(f1, OP1, OP2, f2) \
        TSEQ_FORALL_F_OP_OP_F_impl(f1, OP1, OP2, f2, __restrict__)

#define TSEQ_FORALL_F_OP_OP_F_inplace(f1, OP1, OP2, f2) \
        TSEQ_FORALL_F_OP_OP_F_impl(f1, OP1, OP2, f2,)

// old-style
#define TFOR_ALL_F_OP_OP_F(typeF1, f1, OP1, OP2, typeF2, f2) \
     TSEQ_FORALL_F_OP_OP_F_impl(f1, OP1, OP2, f2, __restrict__)


// ------------
// Member operator : this field f OP s

#define TSEQ_FORALL_F_OP_S_impl(f, OP, s, Unused)                              \
{                                                                              \
    /* Field access */                                                         \
    TFOR_List_ACCESS(f, fP, __restrict__);                                     \
    TFOR_List_LENGTH(f, loop_len);                                             \
                                                                               \
    /* Loop: f OP s */                                                         \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (fP[i]) OP (s);                                                        \
    }                                                                          \
}

#define TSEQ_FORALL_F_OP_S(f, OP, s) \
        TSEQ_FORALL_F_OP_S_impl(f, OP, s, __restrict__)

#define TSEQ_FORALL_F_OP_S_inplace(f, OP, s) \
        TSEQ_FORALL_F_OP_S_impl(f, OP, s,)

// old-style
#define TFOR_ALL_F_OP_S(typeF, f, OP, typeS, s) \
     TSEQ_FORALL_F_OP_S_impl(f, OP, s, __restrict__)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ------------
// Friend operator function : s OP f

#define TSEQ_FORALL_S_OP_F_impl(s, OP, f, Unused)                              \
{                                                                              \
    /* Field access */                                                         \
    TFOR_List_CONST_ACCESS(f, fP, __restrict__);                               \
    TFOR_List_LENGTH(f, loop_len);                                             \
                                                                               \
    /* Loop: s OP f */                                                         \
                                                                               \
    /* pragmas, reduction... */                                                \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (s) OP (fP[i]);                                                        \
    }                                                                          \
}

#define TSEQ_FORALL_S_OP_F(s, OP, f) \
        TSEQ_FORALL_S_OP_F_impl(s, OP, f, __restrict__)

#define TSEQ_FORALL_S_OP_F_inplace(s, OP, f) \
        TSEQ_FORALL_S_OP_F_impl(s, OP, f,)

// old-style
#define TFOR_ALL_S_OP_F(typeS, s, OP, typeF, f) \
     TSEQ_FORALL_S_OP_F_impl(s, OP, f, __restrict__)


// ------------
// Friend operator function : s OP1 f1 OP2 f2

#define TSEQ_FORALL_S_OP_F_OP_F_impl(s, OP1, f1, OP2, f2, Unused)              \
{                                                                              \
    /* Field access */                                                         \
    TFOR_List_ACCESS(f1, f1P, __restrict__);                                   \
    TFOR_List_CONST_ACCESS(f2, f2P, __restrict__);                             \
    TFOR_List_LENGTH(f1, loop_len);                                            \
                                                                               \
    /* Loop: s OP1 f1 OP2 f2 */                                                \
                                                                               \
    /* pragmas, reduction... */                                                \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (s) OP1 (f1P[i]) OP2 (f2P[i]);                                         \
    }                                                                          \
}

#define TSEQ_FORALL_S_OP_F_OP_F(s, OP1, f1, OP2, f2) \
        TSEQ_FORALL_S_OP_F_OP_F_impl(s, OP1, f1, OP2, f2, __restrict__)

#define TSEQ_FORALL_S_OP_F_OP_F_inplace(s, OP1, f1, OP2, f2) \
        TSEQ_FORALL_S_OP_F_OP_F_impl(s, OP1, f1, OP2, f2,)

// old-style
#define TFOR_ALL_S_OP_F_OP_F(typeS, s, OP1, type1, f1, OP2, type2, f2) \
     TSEQ_FORALL_S_OP_F_OP_F_impl(s, OP1, f1, OP2, f2, __restrict__)


// ------------
// Friend operator function : s OP FUNC(f)

#define TSEQ_FORALL_S_OP_FUNC_F_impl(s, OP, FUNC, f, Unused)                   \
{                                                                              \
    /* Field access */                                                         \
    TFOR_List_CONST_ACCESS(f, fP, __restrict__);                               \
    TFOR_List_LENGTH(f, loop_len);                                             \
                                                                               \
    /* Loop: s OP FUNC(f) */                                                   \
                                                                               \
    /* pragmas, reduction... */                                                \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (s) OP FUNC(fP[i]);                                                    \
    }                                                                          \
}

#define TSEQ_FORALL_S_OP_FUNC_F(s, OP, FUNC, f) \
        TSEQ_FORALL_S_OP_FUNC_F_impl(s, OP, FUNC, f, __restrict__)

#define TSEQ_FORALL_S_OP_FUNC_F_inplace(s, OP, FUNC, f) \
        TSEQ_FORALL_S_OP_FUNC_F_impl(s, OP, FUNC, f,)

// old-style
#define TFOR_ALL_S_OP_FUNC_F(typeS, s, OP, FUNC, typeF, f) \
     TSEQ_FORALL_S_OP_FUNC_F_impl(s, OP, FUNC, f, __restrict__)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
