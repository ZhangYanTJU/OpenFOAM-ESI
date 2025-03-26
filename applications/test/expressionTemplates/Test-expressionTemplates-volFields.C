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
            p.expr() + Expression::sqr(p.expr())
        );
        DebugVar(result);
    }
    {
        typedef Expression::GeometricFieldConstRefWrap
        <
            volScalarField
        > volExpr;

        typedef Expression::GeometricFieldConstRefWrap
        <
            surfaceScalarField
        > surfaceExpr;

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
            Expression::lerp<volExpr, surfaceExpr>
            (
                p.expr(),
                mesh.surfaceInterpolation::weights().expr(),
                mesh
            )
        );
        DebugVar(result);
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
