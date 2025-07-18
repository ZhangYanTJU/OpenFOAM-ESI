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
    Test-SphericalTensor2D

Description
    Tests for \c SphericalTensor2D constructors, member functions and operators
    using \c floatScalar, \c doubleScalar, and \c complex base types.

    Cross-checks were obtained from 'NumPy 1.15.1' and 'SciPy 1.1.0' if no
    theoretical cross-check exists (like eigendecomposition relations), and
    were hard-coded for elementwise comparisons.

    For \c complex base type, the cross-checks do only involve zero imag part.

\*---------------------------------------------------------------------------*/

#include "complex.H"
#include "Tensor2D.H"
#include "SymmTensor2D.H"
#include "SphericalTensor2D.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Total number of unit tests
unsigned nTest_ = 0;


// Total number of failed unit tests
unsigned nFail_ = 0;


// Compare two containers elementwise, and print output.
// Do ++nFail_ if values of two objects are not equal within a given tolerance.
// The function is converted from PEP-485.
template<class Type>
void cmp
(
    const word& msg,
    const Type& x,
    const Type& y,
    const scalar relTol = 1e-8,  //<! are values the same within 8 decimals
    const scalar absTol = 0      //<! useful for cmps near zero
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


// Create each constructor of SphericalTensor2D<Type>, and print output
template<class Type>
void test_constructors(Type)
{
    {
        Info<< "# Construct initialized to zero:" << nl;
        const SphericalTensor2D<Type> spT(Zero);
        Info<< spT << endl;
    }
    {
        Info<< "# Construct given VectorSpace of the same rank:" << nl;
        const VectorSpace<SphericalTensor2D<Type>, Type, 1> V(Zero);
        const SphericalTensor2D<Type> spT(V);
        Info<< spT << endl;
    }
    {
        Info<< "# Construct given the component:" << nl;
        const SphericalTensor2D<Type> spT(Type(1));
        Info<< spT << endl;
    }
    {
        Info<< "# Copy construct:" << nl;
        const SphericalTensor2D<Type> spT(Zero);
        const SphericalTensor2D<Type> copyspT(spT);
        Info<< spT << tab << copyspT << endl;
    }
}


// Execute each member function of SphericalTensor2D<Type>, and print output
template<class Type>
void test_member_funcs(Type)
{
    SphericalTensor2D<Type> spT(Type(1));
    const SphericalTensor2D<Type> cspT(Type(-9));

    Info<< "# Operand: " << nl
        << "  SphericalTensor2D = " << spT << endl;


    {
        Info<< "# Component access:" << nl;

        SphericalTensor2D<Type> cpspT(spT.ii());
        cmp("  'SphericalTensor2D' access:", spT, cpspT);

        const SphericalTensor2D<Type> cpcspT(cspT.ii());
        cmp("  'const SphericalTensor2D' access:", cspT, cpcspT);
    }
}


// Execute each global function of SphericalTensor2D<Type>, and print output
template<class Type>
void test_global_funcs(Type)
{
    const SphericalTensor2D<Type> spT(Type(5));

    Info<< "# Operand: " << nl
        << "  SphericalTensor2D = " << spT << endl;


    cmp("  Trace = ", tr(spT), Type(10));
    cmp("  Spherical part = ", sph(spT), spT);
    cmp("  Determinant = ", det(spT), Type(24.99999999999994));
    cmp
    (
        "  Inverse = ",
        inv(spT),
        SphericalTensor2D<Type>(Type(0.2))
    );
}


// Execute each global operator of SphericalTensor2D<Type>, and print output
template<class Type>
void test_global_opers(Type)
{
    const Tensor2D<Type> T
    (
        Type(-1), Type(2),
        Type(4),  Type(5)
    );
    const SymmTensor2D<Type> sT
    (
        Type(-1), Type(2),
                  Type(5)
    );
    const SphericalTensor2D<Type> spT(Type(-2));
    const Vector2D<Type> v(Type(3), Type(2));
    const Type x(4);

    Info<< "# Operands:" << nl
        << "  Tensor2D = " << T << nl
        << "  SymmTensor2D = " << sT << nl
        << "  SphericalTensor2D = " << spT << nl
        << "  Vector2D = " << v << nl
        << "  Type = " << x << endl;


    cmp
    (
        "  Division of Type by SpTensor2D = ",
        (x/spT),
        SphericalTensor2D<Type>(Type(-2))
    );
    cmp
    (
        "  Division of SpTensor2D by Type = ",
        (spT/x),
        SphericalTensor2D<Type>(Type(-0.5))
    );
    cmp
    (
        "  Inner-product of SpTensor2D-SpTensor2D = ",
        (spT & spT),
        SphericalTensor2D<Type>(Type(4))
    );
    cmp
    (
        "  Inner-product of SpTensor2D-Vector2D = ",
        (spT & v),
        Vector2D<Type>(Type(-6), Type(-4)) // Column-vector
    );
    cmp
    (
        "  Inner-product of Vector2D-SpTensor2D = ",
        (v & spT),
        Vector2D<Type>(Type(-6), Type(-4)) // Row-vector
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
        "SphericalTensor2D<float>",
        "SphericalTensor2D<double>",
        "SphericalTensor2D<complex>"
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
