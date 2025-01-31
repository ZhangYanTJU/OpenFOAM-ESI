/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2025 OpenCFD Ltd.
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

Application
    Test-DiagTensor

Description
    Tests for \c DiagTensor constructors, member functions and operators
    using \c floatScalar, \c doubleScalar, and \c complex base types.

    Cross-checks were obtained from 'NumPy 1.15.1' and 'SciPy 1.1.0' if no
    theoretical cross-check exists (like eigendecomposition relations), and
    were hard-coded for elementwise comparisons.

    For \c complex base type, the cross-checks do only involve zero imag part.

\*---------------------------------------------------------------------------*/

#include "complex.H"
#include "Tensor.H"
#include "SymmTensor.H"
#include "SphericalTensor.H"
#include "DiagTensor.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Total number of unit tests
unsigned nTest_ = 0;


// Total number of failed unit tests
unsigned nFail_ = 0;


// Compare two containers elementwise, and print output.
// Do ++nFail_ if two components are not equal within a given tolerance.
// The function is converted from PEP-485
template<class Type>
void cmp
(
    const word& msg,
    const Type& x,
    const Type& y,
    const scalar relTol = 1e-8,
    const scalar absTol = 0
)
{
    const auto notEqual = [=](const auto& a, const auto& b) -> bool
    {
        return
        (
            Foam::max(absTol, relTol*Foam::max(Foam::mag(a), Foam::mag(b)))
          < Foam::mag(a - b)
        );
    };

    unsigned nFail = 0;

    if constexpr (is_vectorspace_v<Type>)
    {
        for (direction i = 0; i < pTraits<Type>::nComponents; ++i)
        {
            if (notEqual(x[i], y[i]))
            {
                ++nFail;
            }
        }
    }
    else
    {
        if (notEqual(x, y))
        {
            ++nFail;
        }
    }

    Info<< msg << x << endl;

    if (nFail)
    {
        Info<< nl
            << "        #### Fail in " << nFail << " comps ####" << nl << endl;
        ++nFail_;
    }
    ++nTest_;
}


// Create each constructor of DiagTensor<Type>, and print output
template<class Type>
void test_constructors(Type)
{
    {
        Info<< "# Construct initialized to zero:" << nl;
        const DiagTensor<Type> dT(Zero);
        Info<< dT << endl;
    }
    {
        Info<< "# Construct given VectorSpace of the same rank:" << nl;
        const VectorSpace<DiagTensor<Type>, Type, 3> V(Zero);
        const DiagTensor<Type> dT(V);
        Info<< dT << endl;
    }
    {
        Info<< "# Construct given the three components:" << nl;
        const DiagTensor<Type> dT
        (
            Type(1),
                    Type(5),
                            Type(-9)
        );
        Info<< dT << endl;
    }
    {
        Info<< "# Copy construct:" << nl;
        const DiagTensor<Type> dT(Zero);
        const DiagTensor<Type> copydT(dT);
        Info<< dT << tab << copydT << endl;
    }
}


// Execute each member function of DiagTensor<Type>, and print output
template<class Type>
void test_member_funcs(Type)
{
    DiagTensor<Type> dT(Type(1), Type(5), Type(-9));
    const DiagTensor<Type> cdT(Type(-9), Type(5), Type(1));

    Info<< "# Operand: " << nl
        << "  DiagTensor = " << dT << endl;


    {
        Info<< "# Component access:" << nl;

        DiagTensor<Type> cpdT(dT.xx(), dT.yy(), dT.zz());
        cmp("  'DiagTensor' access:", dT, cpdT);

        const DiagTensor<Type> cpcdT(cdT.xx(), cdT.yy(), cdT.zz());
        cmp("  'const DiagTensor' access:", cdT, cpcdT);
    }
}


// Execute each global function of DiagTensor<Type>, and print output
template<class Type>
void test_global_funcs(Type)
{
    const Tensor<Type> T
    (
        Type(-1), Type(2), Type(-3),
        Type(4),  Type(5), Type(-6),
        Type(7),  Type(8), Type(-9)
    );
    const SymmTensor<Type> sT
    (
        Type(-1), Type(2), Type(-3),
                  Type(5), Type(-6),
                           Type(-9)
    );
    const DiagTensor<Type> dT(Type(1), Type(5), Type(-9));

    Info<< "# Operands: " << nl
        << "  Tensor = " << T << nl
        << "  SymmTensor = " << sT << nl
        << "  DiagTensor = " << dT << endl;


    cmp("  Trace = ", tr(dT), Type(-3));
    cmp("  Spherical part = ", sph(dT), SphericalTensor<Type>(tr(dT)/Type(3)));
    cmp("  Determinant = ", det(dT), Type(-44.99999999999999));
    cmp
    (
        "  Inverse = ",
        inv(dT),
        DiagTensor<Type>(Type(1), Type(0.2), Type(-0.11111111))
    );
    cmp
    (
        "  Diagonal of Tensor = ",
        diag(T),
        DiagTensor<Type>(Type(-1), Type(5), Type(-9))
    );
    cmp
    (
        "  Diagonal of SymmTensor = ",
        diag(sT),
        DiagTensor<Type>(Type(-1), Type(5), Type(-9))
    );
}


// Execute each global operator of DiagTensor<Type>, and print output
template<class Type>
void test_global_opers(Type)
{
    const Tensor<Type> T
    (
        Type(-1), Type(2), Type(-3),
        Type(4),  Type(5), Type(-6),
        Type(7),  Type(8), Type(-9)
    );
    const SymmTensor<Type> sT
    (
        Type(-1), Type(2), Type(-3),
                  Type(5), Type(-6),
                           Type(-9)
    );
    const DiagTensor<Type> dT(Type(1), Type(5), Type(-9));
    const SphericalTensor<Type> spT(Type(1));
    const Vector<Type> v(Type(3), Type(2), Type(1));
    const Type x(4);

    Info<< "# Operands:" << nl
        << "  Tensor = " << T << nl
        << "  SymmTensor = " << sT << nl
        << "  DiagTensor = " << dT << nl
        << "  SphericalTensor = " << spT << nl
        << "  Vector = " << v << nl
        << "  Type = " << x << endl;


    cmp
    (
        "  Sum of DiagTensor-Tensor = ",
        (dT + T),
        Tensor<Type>
        (
            Type(0), Type(2),  Type(-3),
            Type(4), Type(10), Type(-6),
            Type(7), Type(8),  Type(-18)
        )
    );
    cmp
    (
        "  Sum of Tensor-DiagTensor = ",
        (T + dT),
        Tensor<Type>
        (
            Type(0), Type(2),  Type(-3),
            Type(4), Type(10), Type(-6),
            Type(7), Type(8),  Type(-18)
        )
    );
    cmp
    (
        "  Subtract Tensor from DiagTensor = ",
        (dT - T),
        Tensor<Type>
        (
            Type(2),  Type(-2), Type(3),
            Type(-4), Type(0),  Type(6),
            Type(-7), Type(-8), Type(0)
        )
    );
    cmp
    (
        "  Subtract DiagTensor from Tensor = ",
        (T - dT),
        Tensor<Type>
        (
            Type(-2), Type(2), Type(-3),
            Type(4),  Type(0), Type(-6),
            Type(7),  Type(8), Type(0)
        )
    );
    cmp
    (
        "  Division of Type by DiagTensor = ",
        (x/dT),
        DiagTensor<Type>(Type(4), Type(0.8), Type(-0.44444444))
    );
    cmp
    (
        "  Division of DiagTensor by Type = ",
        (dT/x),
        DiagTensor<Type>(Type(0.25), Type(1.25), Type(-2.25))
    );
    cmp
    (
        "  Division of Vector by DiagTensor = ",
        (v/dT),
        Vector<Type>(Type(3), Type(0.4), Type(-0.11111111))
    );
    cmp
    (
        "  Inner-product of DiagTensor-DiagTensor = ",
        (dT & dT),
        DiagTensor<Type>(Type(1), Type(25), Type(81))
    );
    cmp
    (
        "  Inner-product of DiagTensor-Tensor = ",
        (dT & T),
        Tensor<Type>
        (
            Type(-1),  Type(2),   Type(-3),
            Type(20),  Type(25),  Type(-30),
            Type(-63), Type(-72), Type(81)
        )
    );
    cmp
    (
        "  Inner-product of Tensor-DiagTensor = ",
        (T & dT),
        Tensor<Type>
        (
            Type(-1), Type(10), Type(27),
            Type(4),  Type(25), Type(54),
            Type(7),  Type(40), Type(81)
        )
    );
    cmp
    (
        "  Inner-product of DiagTensor-Vector = ",
        (dT & v),
        Vector<Type>(Type(3), Type(10), Type(-9))
    );
    cmp
    (
        "  Inner-product of Vector-DiagTensor = ",
        (v & dT),
        Vector<Type>(Type(3), Type(10), Type(-9))
    );
}


// Do compile-time recursion over the given types
template<std::size_t I = 0, typename... Tp>
void run_tests(const std::tuple<Tp...>& types, const List<word>& names)
{
    if constexpr (I < sizeof...(Tp))
    {
        const auto& name = names[I];

        Info<< nl << "    ## Test constructors: " << name << " ##" << nl;
        test_constructors(std::get<I>(types));

        Info<< nl << "    ## Test member functions: " << name << " ##" << nl;
        test_member_funcs(std::get<I>(types));

        Info<< nl << "    ## Test global functions: " << name << " ##" << nl;
        test_global_funcs(std::get<I>(types));

        Info<< nl << "    ## Test global operators: " << name << " ##" << nl;
        test_global_opers(std::get<I>(types));

        run_tests<I + 1, Tp...>(types, names);
    }
}


// * * * * * * * * * * * * * * * Main Program  * * * * * * * * * * * * * * * //

int main()
{
    const std::tuple<floatScalar, doubleScalar, complex> types
    (
        std::make_tuple(Zero, Zero, Zero)
    );

    const List<word> typeID
    ({
        "DiagTensor<float>",
        "DiagTensor<double>",
        "DiagTensor<complex>"
    });

    run_tests(types, typeID);


    if (nFail_)
    {
        Info<< nl << "        #### "
            << "Failed in " << nFail_ << " tests "
            << "out of total " << nTest_ << " tests "
            << "####\n" << endl;
        return 1;
    }

    Info<< nl << "        #### Passed all " << nTest_ <<" tests ####\n" << endl;
    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
