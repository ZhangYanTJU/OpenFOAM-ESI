/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2010-2015 Hrvoje Jasak, Wikki Ltd.
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
    Preconditioned Generalised Minimal Residual solver with
    run-time selectable preconditioning

\*---------------------------------------------------------------------------*/

#include "gmresSolver.H"
#include "scalarMatrices.H"
#include "PrecisionAdaptor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gmresSolver, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<gmresSolver>
        addgmresSolverSymMatrixConstructorToTable_;

    lduMatrix::solver::addasymMatrixConstructorToTable<gmresSolver>
        addgmresSolverAsymMatrixConstructorToTable_;

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::gmresSolver::givensRotation
(
    const solveScalar& h,
    const solveScalar& beta,
    solveScalar& c,
    solveScalar& s
) const
{
    if (beta == 0)
    {
        c = 1;
        s = 0;
    }
    else if (mag(beta) > mag(h))
    {
        scalar tau = -h/beta;
        s = 1.0/Foam::sqrt(1.0 + sqr(tau));
        c = s*tau;
    }
    else
    {
        scalar tau = -beta/h;
        c = 1.0/Foam::sqrt(1.0 + sqr(tau));
        s = c*tau;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix and solver data stream
Foam::gmresSolver::gmresSolver
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& dict
)
:
    lduMatrix::solver
    (
        fieldName,
        matrix,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaces,
        dict
    ),
    preconPtr_
    (
        lduMatrix::preconditioner::New
        (
            *this,
            dict
        )
    ),
    nDirs_(dict.getLabel("nDirections"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::solverPerformance Foam::gmresSolver::scalarSolve
(
    solveScalarField& x,
    const solveScalarField& b,
    const direction cmpt
) const
{
    // Initialise the solverPerformance object to track solver performance
    solverPerformance solverPerf
    (
        preconPtr_->type() + typeName, fieldName()
    );

    solveScalarField wA(x.size());
    solveScalarField rA(x.size());

    // Calculate initial residual
    matrix_.Amul(wA, x, interfaceBouCoeffs_, interfaces_, cmpt);

    // Use rA as scratch space when calculating the normalisation factor
    solveScalar normFactor = this->normFactor(x, b, wA, rA);

    if (lduMatrix::debug >= 2)
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }

    // Calculate residual
    forAll (rA, i)
    {
        rA[i] = b[i] - wA[i];
    }

    solverPerf.initialResidual() =
        gSumMag(rA, matrix().mesh().comm())
       /normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    // Note: GMRES cannot be forced to do minIter sweeps
    // if the residual is zero, due to algorithmic reasons
    // HJ, 22/Aug/2012
    if (!solverPerf.checkConvergence(tolerance_, relTol_, log_))
    {
        typedef SquareMatrix<solveScalar> solveScalarSquareMatrix;

        // Create the Hesenberg matrix
        solveScalarSquareMatrix H(nDirs_, 0);

        // Create y and b for Hessenberg matrix
        solveScalarField yh(nDirs_, 0);
        solveScalarField bh(nDirs_ + 1, 0);

        // Givens rotation vectors
        solveScalarField c(nDirs_, 0);
        solveScalarField s(nDirs_, 0);

        // Allocate Krylov space vectors
        FieldField<Field, solveScalar> V(nDirs_ + 1);

        forAll (V, i)
        {
            V.set(i, new solveScalarField(x.size(), 0));
        }

        do
        {
            // Execute preconditioning
            preconPtr_->precondition(wA, rA, cmpt);

            // Calculate beta and scale first vector
            solveScalar beta = Foam::sqrt(gSumSqr(wA, matrix().mesh().comm()));

            // Set initial rhs and bh[0] = beta
            bh = 0;
            bh[0] = beta;

            for (label i = 0; i < nDirs_; ++i)
            {
                // Set search direction
                V[i] = wA;
                V[i] /= beta;

                // Arnoldi's method
                matrix_.Amul(rA, V[i], interfaceBouCoeffs_, interfaces_, cmpt);

                // Execute preconditioning
                preconPtr_->precondition(wA, rA, cmpt);

                for (label j = 0; j <= i; ++j)
                {
                    beta = gSumProd(wA, V[j], matrix().mesh().comm());

                    H[j][i] = beta;

                    forAll (wA, wI)
                    {
                        wA[wI] -= beta*V[j][wI];
                    }
                }

                beta = Foam::sqrt(gSumSqr(wA, matrix().mesh().comm()));

                // Apply previous Givens rotations to new column of H.
                for (label j = 0; j < i; ++j)
                {
                    const solveScalar Hji = H[j][i];
                    H[j][i] = c[j]*Hji - s[j]*H[j + 1][i];
                    H[j + 1][i] = s[j]*Hji + c[j]*H[j + 1][i];
                }

                // Apply Givens rotation to current row.
                givensRotation(H[i][i], beta, c[i], s[i]);

                const solveScalar bhi = bh[i];
                bh[i] = c[i]*bhi - s[i]*bh[i + 1];
                bh[i + 1] = s[i]*bhi + c[i]*bh[i + 1];
                H[i][i] = c[i]*H[i][i] - s[i]*beta;
            }

            // Back substitute to solve Hy = b
            for (label i = nDirs_ - 1; i >= 0; i--)
            {
                solveScalar sum = bh[i];

                for (label j = i + 1; j < nDirs_; ++j)
                {
                    sum -= H[i][j]*yh[j];
                }

                yh[i] = sum/H[i][i];
            }

            // Update solution
            for (label i = 0; i < nDirs_; ++i)
            {
                const solveScalarField& Vi = V[i];
                const solveScalar& yi = yh[i];

                forAll (x, psiI)
                {
                    x[psiI] += yi*Vi[psiI];
                }
            }

            // Re-calculate the residual
            matrix_.Amul(wA, x, interfaceBouCoeffs_, interfaces_, cmpt);

            forAll (rA, raI)
            {
                rA[raI] = b[raI] - wA[raI];
            }

            solverPerf.finalResidual() =
                gSumMag(rA, matrix().mesh().comm())
               /normFactor;
            solverPerf.nIterations()++;
        } while
        (
            (
              ++solverPerf.nIterations() < maxIter_
            && !solverPerf.checkConvergence(tolerance_, relTol_, log_)
            )
         || solverPerf.nIterations() < minIter_
        );
    }


    if (preconPtr_)
    {
        preconPtr_->setFinished(solverPerf);
    }

    matrix().setResidualField
    (
        ConstPrecisionAdaptor<scalar, solveScalar>(rA)(),
        fieldName_,
        false
    );

    return solverPerf;
}


Foam::solverPerformance Foam::gmresSolver::solve
(
    scalarField& psi_s,
    const scalarField& source,
    const direction cmpt
) const
{
    PrecisionAdaptor<solveScalar, scalar> tpsi(psi_s);
    return scalarSolve
    (
        tpsi.ref(),
        ConstPrecisionAdaptor<solveScalar, scalar>(source)(),
        cmpt
    );
}


// ************************************************************************* //
