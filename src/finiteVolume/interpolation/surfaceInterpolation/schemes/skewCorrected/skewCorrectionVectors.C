/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "skewCorrectionVectors.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(skewCorrectionVectors, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::skewCorrectionVectors::skewCorrectionVectors(const fvMesh& mesh)
:
    MeshObject_type(mesh),
    skew_(false),
    skewCorrectionVectors_
    (
        IOobject
        (
            "skewCorrectionVectors",
            mesh_.pointsInstance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimless
    )
{
    calcSkewCorrectionVectors();
}


Foam::skewCorrectionVectors::~skewCorrectionVectors()
{}


void Foam::skewCorrectionVectors::calcSkewCorrectionVectors()
{
    DebugInFunction << "Calculating skew correction vectors" << nl;

    // Set local references to mesh data
    const volVectorField& C = mesh_.C();
    const surfaceVectorField& Cf = mesh_.Cf();
    const surfaceVectorField& Sf = mesh_.Sf();

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        vector d = C[nei] - C[own];
        vector Cpf = Cf[facei] - C[own];

        skewCorrectionVectors_[facei] =
            Cpf - ((Sf[facei] & Cpf)/(Sf[facei] & d))*d;
    }

    surfaceVectorField::Boundary& skewCorrVecsBf =
        skewCorrectionVectors_.boundaryFieldRef();

    forAll(skewCorrVecsBf, patchi)
    {
        fvsPatchVectorField& patchSkewCorrVecs = skewCorrVecsBf[patchi];

        if (!patchSkewCorrVecs.coupled())
        {
            patchSkewCorrVecs = Zero;
        }
        else
        {
            const fvPatch& p = patchSkewCorrVecs.patch();
            const labelUList& faceCells = p.faceCells();
            const vectorField& patchFaceCentres = Cf.boundaryField()[patchi];
            const vectorField& patchSf = Sf.boundaryField()[patchi];
            const vectorField patchD(p.delta());

            forAll(p, patchFacei)
            {
                vector Cpf =
                    patchFaceCentres[patchFacei] - C[faceCells[patchFacei]];

                patchSkewCorrVecs[patchFacei] =
                    Cpf
                  - (
                        (patchSf[patchFacei] & Cpf)/
                        (patchSf[patchFacei] & patchD[patchFacei])
                    )*patchD[patchFacei];
            }
        }
    }

    scalar skewCoeff = 0.0;

    if (Sf.primitiveField().size())
    {
        skewCoeff =
            max(mag(skewCorrectionVectors_)*mesh_.deltaCoeffs()).value();
    }

    DebugInFunction << "skew coefficient = " << skewCoeff << nl;

    if (skewCoeff < 1e-5)
    {
        skew_ = false;
    }
    else
    {
        skew_ = true;
    }

    if (debug)
    {
        Info<< "    Finished constructing skew correction vectors" << endl;
    }
}


bool Foam::skewCorrectionVectors::movePoints()
{
    calcSkewCorrectionVectors();
    return true;
}


// ************************************************************************* //
