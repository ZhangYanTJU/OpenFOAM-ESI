/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    volPointInterpolationTest

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "argList.H"
#include "fvMesh.H"
#include "fvCFD.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volScalarField> someFunction(const volScalarField& fld)
{
    return fld*1.0001;
}


template<class Type>
void fusedGaussFvmLaplacian
(
    fvMatrix<Type>& fvm,
    const surfaceInterpolationScheme<scalar>& interpGammaScheme,
    const fv::snGradScheme<Type>& snGradScheme,
    const GeometricField<scalar, fvPatchField, volMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    // Replacement for gaussLaplacianScheme::fvmLaplacian with scalar gamma
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> surfaceType;

    const auto& mesh = vf.mesh();

    // Expression for weights
    const auto weights = interpGammaScheme.weights(gamma).expr();

    // Expression for gamma_face * magSf
    const auto gammaMagSf =
        Expression::lerp(gamma.expr(), weights, mesh)
      * mesh.magSf().expr();

    // Expression for deltaCoeffs
    const auto deltaCoeffs = snGradScheme.deltaCoeffs(vf).expr();

    // Construct matrix
    Expression::fvmLaplacianUncorrected(fvm, gammaMagSf, deltaCoeffs);

    if (snGradScheme.corrected())
    {
        // Wrap correction
        const auto corr(snGradScheme.correction(vf).expr());
        const auto V = mesh.V().expr();

        if (mesh.fluxRequired(vf.name()))
        {
            fvm.faceFluxCorrectionPtr() = std::make_unique<surfaceType>
            (
                IOobject
                (
                    "faceFluxCorr",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh,
                gamma.dimensions()
               *mesh.magSf().dimensions()
               *corr.data().dimensions()
            );
            auto& faceFluxCorr = *fvm.faceFluxCorrectionPtr();
            faceFluxCorr = gammaMagSf*corr;

            fvm.source() =
                fvm.source().expr()
              - (
                    V * fvc::div
                    (
                        faceFluxCorr
                    )().primitiveField().expr()
                );
        }
        else
        {
            // Temporary field
            surfaceType faceFluxCorr
            (
                IOobject
                (
                    "faceFluxCorr",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh,
                gamma.dimensions()
               *mesh.magSf().dimensions()
               *corr.data().dimensions()
            );
            faceFluxCorr = gammaMagSf*corr;

            fvm.source() =
                fvm.source().expr()
              - (
                    V * fvc::div
                    (
                        faceFluxCorr
                    )().primitiveField().expr()
                );
        }
    }
}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh.thisDb(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Expresions of volFields
    {
        volScalarField result
        (
            IOobject
            (
                "result",
                runTime.timeName(),
                mesh.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensionedScalar(p.dimensions(), 0)
        );
        auto expression = someFunction(p).expr() + someFunction(p).expr();
        result = expression;
        DebugVar(result);
    }

    // Expresions of volFields
    {
        volScalarField result
        (
            "result",
            mesh,
            sqr(p.expr() + p.expr())
        );
        DebugVar(result);
    }
    {
        // Fill p with some values
        forAll(p, celli)
        {
            p[celli] = celli;
        }
        p.correctBoundaryConditions();

        // Interpolate to surface field
        surfaceScalarField result
        (
            "result",
            mesh,
            Expression::lerp    //<volExpr, surfaceExpr>
            (
                p.expr(),
                mesh.surfaceInterpolation::weights().expr(),
                mesh
            )
        );
        DebugVar(result);
    }
    {
        // For testing as a replacement of laplacian weights
        const volScalarField gamma
        (
            IOobject
            (
                "gamma",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensionedScalar(dimless, 1.0)
        );

        fvMatrix<scalar> fvm
        (
            p,
            gamma.dimensions()*mesh.magSf().dimensions()*p.dimensions()
        );

        const linear<scalar> interpGammaScheme(mesh);
        const fv::correctedSnGrad<scalar> snGradScheme(mesh);

        fusedGaussFvmLaplacian
        (
            fvm,
            interpGammaScheme,
            snGradScheme,
            gamma,
            p
        );

        DebugVar(fvm.source());
    }



    // Expressions of fvMatrix
    {
        tmp<fvMatrix<scalar>> tm0(fvm::laplacian(p));
        const fvMatrix<scalar>& m0 = tm0();
        DebugVar(m0.dimensions());

        tmp<fvMatrix<scalar>> tm1(fvm::laplacian(p));
        const fvMatrix<scalar>& m1 = tm1();
        DebugVar(m1.dimensions());

        fvMatrix<scalar> m2(p, m0.expr() + m1.expr());
        DebugVar(m2.dimensions());
    }

    return 0;
}


// ************************************************************************* //
