/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "SVD.H"
#include "scalarList.H"
#include "scalarMatrices.H"
#include "ListOps.H"

#include "tensor.H"
#include "diagTensor.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Same as implementation in scalarMatrices, but for tensor/diagTensor/tensor
static tensor do_multiply
(
    const tensor& A,
    const diagTensor& B,
    const tensor& C
)
{
    constexpr direction size = 3;

    tensor ans(Foam::zero{});

    for (direction i = 0; i < size; ++i)
    {
        for (direction g = 0; g < size; ++g)
        {
            for (direction l = 0; l < size; ++l)
            {
                ans(i, g) += C(l, g)*A(i, l)*B[l];
            }
        }
    }

    return ans;
}


template<class MatrixType, class DiagMatrixType, class WorkArrayType>
static std::pair<label, bool> SVDcomp
(
    // input
    const MatrixType& A,

    // input
    const scalar minCondition,

    // output
    MatrixType& U_,

    // output
    MatrixType& V_,

    // output
    DiagMatrixType& S_,

    // scratch space
    WorkArrayType& rv1
)
{
    label nZeros_(0);
    bool converged_(true);

    // SVDcomp to find U_, V_ and S_ - the singular values

    U_ = A;

    const label Un = U_.n();
    const label Um = U_.m();

    scalar g = 0;
    scalar scale = 0;
    scalar s = 0;
    scalar anorm = 0;
    label l = 0;


    const auto sign = [](scalar a, scalar b) -> scalar
    {
        //return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
        return b >= 0 ? a : -a;
    };


    for (label i=0; i<Un; i++)
    {
        l = i + 2;
        rv1[i] = scale*g;
        g = s = scale = 0;

        if (i < Um)
        {
            for (label k=i; k<Um; k++)
            {
                scale += mag(U_(k, i));
            }

            if (scale != 0)
            {
                for (label k=i; k<Um; k++)
                {
                    U_(k, i) /= scale;
                    s += U_(k, i)*U_(k, i);
                }

                scalar f = U_(i, i);
                g = -sign(sqrt(s), f);
                scalar h = f*g - s;
                U_(i, i) = f - g;

                for (label j=l-1; j<Un; j++)
                {
                    s = 0;
                    for (label k=i; k<Um; k++)
                    {
                        s += U_(k, i)*U_(k, j);
                    }

                    f = s/h;
                    for (label k=i; k<Um; k++)
                    {
                        U_(k, j) += f*U_(k, i);
                    }
                }

                for (label k=i; k<Um; k++)
                {
                    U_(k, i) *= scale;
                }
            }
        }

        S_[i] = scale*g;

        g = s = scale = 0;

        if (i+1 <= Um && i != Un)
        {
            for (label k=l-1; k<Un; k++)
            {
                scale += mag(U_(i, k));
            }

            if (scale != 0)
            {
                for (label k=l-1; k<Un; k++)
                {
                    U_(i, k) /= scale;
                    s += U_(i, k)*U_(i, k);
                }

                scalar f = U_(i, l-1);
                g = -sign(sqrt(s), f);
                scalar h = f*g - s;
                U_(i, l-1) = f - g;

                for (label k=l-1; k<Un; k++)
                {
                    rv1[k] = U_(i, k)/h;
                }

                for (label j=l-1; j<Um; j++)
                {
                    s = 0;
                    for (label k=l-1; k<Un; k++)
                    {
                        s += U_(j, k)*U_(i, k);
                    }

                    for (label k=l-1; k<Un; k++)
                    {
                        U_(j, k) += s*rv1[k];
                    }
                }
                for (label k=l-1; k<Un; k++)
                {
                    U_(i, k) *= scale;
                }
            }
        }

        anorm = max(anorm, mag(S_[i]) + mag(rv1[i]));
    }

    anorm *= SMALL;

    for (label i=Un-1; i >= 0; i--)
    {
        if (i < Un-1)
        {
            if (g*U_(i, l) != 0)
            {
                for (label j=l; j<Un; j++)
                {
                    V_(j, i) = (U_(i, j)/U_(i, l))/g;
                }

                for (label j=l; j<Un; j++)
                {
                    s = 0;
                    for (label k=l; k<Un; k++)
                    {
                        s += U_(i, k)*V_(k, j);
                    }

                    for (label k=l; k<Un; k++)
                    {
                        V_(k, j) += s*V_(k, i);
                    }
                }
            }

            for (label j=l; j<Un;j++)
            {
                V_(i, j) = V_(j, i) = 0;
            }
        }

        V_(i, i) = 1;
        g = rv1[i];
        l = i;
    }

    for (label i=min(Un, Um) - 1; i>=0; i--)
    {
        l = i+1;
        g = S_[i];

        for (label j=l; j<Un; j++)
        {
            U_(i, j) = 0;
        }

        if (g*U_(i, i) != 0)
        {
            g = 1.0/g;

            for (label j=l; j<Un; j++)
            {
                s = 0;
                for (label k=l; k<Um; k++)
                {
                    s += U_(k, i)*U_(k, j);
                }

                scalar f = (s/U_(i, i))*g;

                for (label k=i; k<Um; k++)
                {
                    U_(k, j) += f*U_(k, i);
                }
            }

            for (label j=i; j<Um; j++)
            {
                U_(j, i) *= g;
            }
        }
        else
        {
            for (label j=i; j<Um; j++)
            {
                U_(j, i) = 0;
            }
        }

        ++U_(i, i);
    }

    for (label k=Un-1; k >= 0; k--)
    {
        for (label its = 0; its < 30; its++)
        {
            bool flag = true;

            label mn;
            for (l = k; l >= 0; l--)
            {
                mn = l - 1;

                if (l == 0 || mag(rv1[l]) <= anorm)
                {
                    flag = false;
                    break;
                }

                if (mag(S_[mn]) <= anorm)
                {
                    break;
                }
            }

            if (flag)
            {
                scalar c = 0;
                s = 1;
                for (label i=l; i<k+1; i++)
                {
                    scalar f = s*rv1[i];
                    rv1[i] = c*rv1[i];

                    if (mag(f) <= anorm)
                    {
                        break;
                    }

                    g = S_[i];
                    scalar h = sqrtSumSqr(f, g);
                    S_[i] = h;
                    h = 1.0/h;
                    c = g*h;
                    s = -f*h;

                    for (label j=0; j<Um; j++)
                    {
                        scalar y = U_(j, mn);
                        scalar z = U_(j, i);
                        U_(j, mn) = y*c + z*s;
                        U_(j, i) = z*c - y*s;
                    }
                }
            }

            scalar z = S_[k];

            if (l == k)
            {
                if (z < 0)
                {
                    S_[k] = -z;
                    for (label j=0; j<Un; j++)
                    {
                        V_(j, k) = -V_(j, k);
                    }
                }
                break;
            }

            if (its == 29)
            {
                converged_ = false;
            }

            scalar x = S_[l];
            mn = k-1;
            scalar y = S_[mn];
            g = rv1[mn];
            scalar h = rv1[k];
            scalar f = ((y - z)*(y + z) + (g - h)*(g + h))/(2*h*y);
            g = sqrtSumSqr(f, scalar(1));
            f = ((x - z)*(x + z) + h*((y/(f + sign(g, f))) - h))/x;
            scalar c = 1;
            s = 1;

            for (label j=l; j <= mn; j++)
            {
                label i = j + 1;
                g = rv1[i];
                y = S_[i];
                h = s*g;
                g = c*g;
                scalar z = sqrtSumSqr(f, h);
                rv1[j] = z;
                c = f/z;
                s = h/z;
                f = x*c + g*s;
                g = g*c - x*s;
                h = y*s;
                y *= c;

                for (label jj = 0; jj < Un; jj++)
                {
                    x = V_(jj, j);
                    z = V_(jj, i);
                    V_(jj, j) = x*c + z*s;
                    V_(jj, i) = z*c - x*s;
                }

                z = sqrtSumSqr(f, h);
                S_[j] = z;
                if (z)
                {
                    z = 1.0/z;
                    c = f*z;
                    s = h*z;
                }
                f = c*g + s*y;
                x = c*y - s*g;

                for (label jj=0; jj < Um; jj++)
                {
                    y = U_(jj, j);
                    z = U_(jj, i);
                    U_(jj, j) = y*c + z*s;
                    U_(jj, i) = z*c - y*s;
                }
            }
            rv1[l] = 0;
            rv1[k] = f;
            S_[k] = x;
        }
    }

    // Zero singular values that are less than minCondition*maxS
    const scalar minS = minCondition*S_[findMax(S_)];
    for (auto& val : S_)
    {
        if (val <= minS)
        {
            val = 0;
            ++nZeros_;
        }
    }

    return { nZeros_, converged_ };
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructor  * * * * * * * * * * * * * * //

Foam::SVD::SVD(const scalarRectangularMatrix& A, const scalar minCondition)
:
    U_(),
    V_(A.n(), A.n()),
    S_(A.n()),
    converged_(true),
    nZeros_(0)
{
    scalarList rv1(A.n());

    // SVDcomp to find U_, V_ and S_ - the singular values

    auto status = SVDcomp(A, minCondition, U_, V_, S_, rv1);

    nZeros_ = status.first;
    converged_ = status.second;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::SVD::minNonZeroS() const
{
    scalar minS = S_[0];
    for (label i=1; i<S_.size(); i++)
    {
        scalar s = S_[i];
        if (s > VSMALL && s < minS) minS = s;
    }
    return minS;
}


Foam::scalarRectangularMatrix Foam::SVD::VSinvUt() const
{
    scalarRectangularMatrix VSinvUt;
    multiply(VSinvUt, V_, Foam::inv(S_), U_.T());
    return VSinvUt;
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::scalarRectangularMatrix Foam::SVD::pinv
(
    const scalarRectangularMatrix& A,
    const scalar minCondition
)
{
    SVD svd(A, minCondition);
    return svd.VSinvUt();
}


Foam::Tensor<Foam::scalar>
Foam::SVD::pinv(const Tensor<scalar>& A, const scalar minCondition)
{
    // SVDcomp to find U_, V_ and S_ - the singular values

    tensor U_;
    tensor V_;
    diagTensor S_;
    FixedList<scalar, 3> rv1;

    SVDcomp(A, minCondition, U_, V_, S_, rv1);

    return do_multiply(V_, Foam::inv(S_), U_.T());
}


// ************************************************************************* //
