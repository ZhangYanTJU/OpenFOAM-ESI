/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "areaFields.H"
#include "edgeFields.H"
#include "facGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

inline tmp<areaScalarField> limiter(const areaScalarField& phi)
{
    return phi;
}

inline tmp<areaScalarField> limiter(const areaVectorField& phi)
{
    return magSqr(phi);
}

inline tmp<areaScalarField> limiter(const areaTensorField& phi)
{
    return magSqr(phi);
}

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class NVDweight>
Foam::tmp<Foam::edgeScalarField> Foam::faNVDscheme<Type,NVDweight>::weights
(
    const GeometricField<Type, faPatchField, areaMesh>& phi
) const
{
    const faMesh& mesh = this->mesh();

    auto tWeightingFactors = tmp<edgeScalarField>::New
    (
        mesh.edgeInterpolation::weights()
    );
    edgeScalarField& weightingFactors = tWeightingFactors.ref();

    scalarField& weights = weightingFactors.primitiveFieldRef();

    tmp<areaScalarField> tvf = limiter(phi);
    const areaScalarField& vf = tvf();

    tmp<areaVectorField> tgradc = fac::grad(vf);
    const areaVectorField& gradc = tgradc.cref();

//     edgeVectorField d
//     (
//         mesh.Le()
//        /(mesh.magLe()*mesh.edgeInterpolation::deltaCoeffs())
//     );

//     if (!mesh.orthogonal())
//     {
//         d -=
//             mesh.edgeInterpolation::nonOrthCorrectionVector()
//            /mesh.edgeInterpolation::deltaCoeffs();
//     }

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const vectorField& n = mesh.faceAreaNormals().internalField();
    const vectorField& c = mesh.areaCentres().internalField();

    forAll(weights, edge)
    {
        vector d(c[neighbour[edge]] - c[owner[edge]]);

        if (edgeFlux_[edge] > 0)
        {
            d.removeCollinear(n[owner[edge]]);
        }
        else
        {
            d.removeCollinear(n[neighbour[edge]]);
        }
        d.normalise();
        d *= mesh.edgeInterpolation::lPN().internalField()[edge];

        // Do not allow any mag(val) < SMALL
        if (mag(d) < SMALL)
        {
            d = vector::uniform(SMALL);
        }

        weights[edge] =
            this->weight
            (
                weights[edge],
                edgeFlux_[edge],
                vf[owner[edge]],
                vf[neighbour[edge]],
                gradc[owner[edge]],
                gradc[neighbour[edge]],
                d
            );
    }


    auto& bWeights = weightingFactors.boundaryFieldRef();

    forAll(bWeights, patchi)
    {
        if (bWeights[patchi].coupled())
        {
            scalarField& pWeights = bWeights[patchi];

            const scalarField& pEdgeFlux = edgeFlux_.boundaryField()[patchi];

            scalarField pVfP(vf.boundaryField()[patchi].patchInternalField());

            scalarField pVfN(vf.boundaryField()[patchi].patchNeighbourField());

            vectorField pGradcP
            (
                gradc.boundaryField()[patchi].patchInternalField()
            );

            vectorField pGradcN
            (
                gradc.boundaryField()[patchi].patchNeighbourField()
            );

            vectorField CP
            (
                mesh.areaCentres().boundaryField()[patchi].patchInternalField()
            );

            vectorField CN
            (
                mesh.areaCentres().boundaryField()[patchi]
                .patchNeighbourField()
            );

            vectorField nP
            (
                mesh.faceAreaNormals().boundaryField()[patchi]
               .patchInternalField()
            );

            vectorField nN
            (
                mesh.faceAreaNormals().boundaryField()[patchi]
               .patchNeighbourField()
            );

            scalarField pLPN
            (
                mesh.edgeInterpolation::lPN().boundaryField()[patchi]
            );

            forAll(pWeights, edgei)
            {
                vector d(CN[edgei] - CP[edgei]);

                if (pEdgeFlux[edgei] > 0)
                {
                    d.removeCollinear(nP[edgei]);
                }
                else
                {
                    d.removeCollinear(nN[edgei]);
                }
                d.normalise();
                d *= pLPN[edgei];

                // Do not allow any mag(val) < SMALL
                if (mag(d) < SMALL)
                {
                    d = vector::uniform(SMALL);
                }

                pWeights[edgei] =
                    this->weight
                    (
                        pWeights[edgei],
                        pEdgeFlux[edgei],
                        pVfP[edgei],
                        pVfN[edgei],
                        pGradcP[edgei],
                        pGradcN[edgei],
                        d
                    );
            }
        }
    }

    return tWeightingFactors;
}


// ************************************************************************* //
