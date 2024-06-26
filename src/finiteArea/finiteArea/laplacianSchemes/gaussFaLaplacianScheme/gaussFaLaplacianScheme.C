/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

#include "gaussFaLaplacianScheme.H"
#include "facDiv.H"
#include "faMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<faMatrix<Type>>
gaussLaplacianScheme<Type>::famLaplacian
(
    const edgeScalarField& gamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<edgeScalarField> tdeltaCoeffs = this->tlnGradScheme_().deltaCoeffs(vf);
    const edgeScalarField& deltaCoeffs = tdeltaCoeffs();

    const edgeScalarField gammaMagSf(gamma*this->mesh().magLe());

    tmp<faMatrix<Type>> tfam
    (
        new faMatrix<Type>
        (
            vf,
            deltaCoeffs.dimensions()*gammaMagSf.dimensions()*vf.dimensions()
        )
    );
    faMatrix<Type>& fam = tfam.ref();

    fam.upper() = deltaCoeffs.internalField()*gammaMagSf.internalField();
    fam.negSumDiag();

    forAll(fam.psi().boundaryField(), patchI)
    {
        const faPatchField<Type>& psf = fam.psi().boundaryField()[patchI];
        const faePatchScalarField& patchGamma =
            gammaMagSf.boundaryField()[patchI];

        fam.internalCoeffs()[patchI] = patchGamma*psf.gradientInternalCoeffs();
        fam.boundaryCoeffs()[patchI] =
            -patchGamma*psf.gradientBoundaryCoeffs();
    }

    if (this->tlnGradScheme_().corrected())
    {
        if (this->mesh().fluxRequired(vf.name()))
        {
            fam.faceFluxCorrectionPtr() = std::make_unique
            <
                GeometricField<Type, faePatchField, edgeMesh>
            >
            (
                gammaMagSf*this->tlnGradScheme_().correction(vf)
            );

            fam.source() -=
                this->mesh().S()
               *fac::div
                (
                    *fam.faceFluxCorrectionPtr()
                )().internalField();
        }
        else
        {
            fam.source() -=
                this->mesh().S()
               *fac::div
                (
                    gammaMagSf*this->tlnGradScheme_().correction(vf)
                )().internalField();
        }
    }

    return tfam;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
gaussLaplacianScheme<Type>::facLaplacian
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> tLaplacian
    (
        fac::div(this->tlnGradScheme_().lnGrad(vf)*vf.mesh().magLe())
    );

    tLaplacian.ref().rename("laplacian(" + vf.name() + ')');

    return tLaplacian;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
gaussLaplacianScheme<Type>::facLaplacian
(
    const edgeScalarField& gamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> tLaplacian
    (
        fac::div(gamma*this->tlnGradScheme_().lnGrad(vf)*vf.mesh().magLe())
    );

    tLaplacian.ref().rename
    (
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'
    );

    return tLaplacian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
