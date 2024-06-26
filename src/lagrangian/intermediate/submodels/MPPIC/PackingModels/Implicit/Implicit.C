/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "Implicit.H"
#include "fixedValueFvsPatchField.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvcReconstruct.H"
#include "volPointInterpolation.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PackingModels::Implicit<CloudType>::Implicit
(
    const dictionary& dict,
    CloudType& owner
)
:
    PackingModel<CloudType>(dict, owner, typeName),
    alpha_
    (
        IOobject
        (
            IOobject::scopedName(this->owner().name(), "alpha"),
            this->owner().db().time().timeName(),
            this->owner().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->owner().mesh(),
        dimensionedScalar(dimless, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    phiCorrect_(nullptr),
    uCorrect_(nullptr),
    applyLimiting_(this->coeffDict().lookup("applyLimiting")),
    applyGravity_(this->coeffDict().lookup("applyGravity")),
    alphaMin_(this->coeffDict().getScalar("alphaMin")),
    rhoMin_(this->coeffDict().getScalar("rhoMin"))
{
    alpha_ = this->owner().theta();
    alpha_.oldTime();
}


template<class CloudType>
Foam::PackingModels::Implicit<CloudType>::Implicit
(
    const Implicit<CloudType>& cm
)
:
    PackingModel<CloudType>(cm),
    alpha_(cm.alpha_),
    phiCorrect_(cm.phiCorrect_()),
    uCorrect_(cm.uCorrect_()),
    applyLimiting_(cm.applyLimiting_),
    applyGravity_(cm.applyGravity_),
    alphaMin_(cm.alphaMin_),
    rhoMin_(cm.rhoMin_)
{
    alpha_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PackingModels::Implicit<CloudType>::~Implicit()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::PackingModels::Implicit<CloudType>::cacheFields(const bool store)
{
    PackingModel<CloudType>::cacheFields(store);

    if (store)
    {
        const fvMesh& mesh = this->owner().mesh();
        const dimensionedScalar deltaT = this->owner().db().time().deltaT();
        const word& cloudName = this->owner().name();

        const dimensionedVector& g = this->owner().g();
        const volScalarField& rhoc = this->owner().rho();

        const auto& rhoAverage =
            mesh.lookupObject<AveragingMethod<scalar>>
            (
                IOobject::scopedName(cloudName, "rhoAverage")
            );
        const auto& uAverage =
            mesh.lookupObject<AveragingMethod<vector>>
            (
                IOobject::scopedName(cloudName, "uAverage")
            );
        const auto& uSqrAverage =
            mesh.lookupObject<AveragingMethod<scalar>>
            (
                IOobject::scopedName(cloudName, "uSqrAverage")
            );

        mesh.setFluxRequired(alpha_.name());

        // Property fields
        // ~~~~~~~~~~~~~~~

        // volume fraction field
        alpha_ = max(this->owner().theta(), alphaMin_);
        alpha_.correctBoundaryConditions();

        // average density
        auto trho = volScalarField::New
        (
            IOobject::scopedName(cloudName, "rho"),
            mesh,
            dimensionedScalar(dimDensity, Zero),
            fvPatchFieldBase::zeroGradientType()
        );
        auto& rho = trho.ref();

        rho.primitiveFieldRef() = max(rhoAverage.primitiveField(), rhoMin_);
        rho.correctBoundaryConditions();

        // Stress field
        // ~~~~~~~~~~~~

        // stress derivative wrt volume fraction
        auto ttauPrime = volScalarField::New
        (
            IOobject::scopedName(cloudName, "tauPrime"),
            mesh,
            dimensionedScalar(dimPressure, Zero),
            fvPatchFieldBase::zeroGradientType()
        );
        auto& tauPrime = ttauPrime.ref();

        tauPrime.primitiveFieldRef() =
            this->particleStressModel_->dTaudTheta
            (
                alpha_.primitiveField(),
                rho.primitiveField(),
                uSqrAverage.primitiveField()
            )();

        tauPrime.correctBoundaryConditions();


        // Gravity flux
        // ~~~~~~~~~~~~

        tmp<surfaceScalarField> phiGByA;

        if (applyGravity_)
        {
            phiGByA = surfaceScalarField::New
            (
                "phiGByA",
                IOobject::NO_REGISTER,
                deltaT*(g & mesh.Sf())*fvc::interpolate(1.0 - rhoc/rho)
            );
        }


        // Implicit solution for the volume fraction
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        surfaceScalarField tauPrimeByRhoAf
        (
            "tauPrimeByRhoAf",
            fvc::interpolate(deltaT*tauPrime/rho)
        );

        fvScalarMatrix alphaEqn
        (
            fvm::ddt(alpha_)
          - fvc::ddt(alpha_)
          - fvm::laplacian(tauPrimeByRhoAf, alpha_)
        );

        if (applyGravity_)
        {
            alphaEqn += fvm::div(phiGByA(), alpha_);
        }

        alphaEqn.solve();


        // Generate correction fields
        // ~~~~~~~~~~~~~~~~~

        // correction volumetric flux
        phiCorrect_ = surfaceScalarField::New
        (
            IOobject::scopedName(cloudName, "phiCorrect"),
            IOobject::NO_REGISTER,
            alphaEqn.flux()/fvc::interpolate(alpha_)
        );

        // limit the correction flux
        if (applyLimiting_)
        {
            auto tU = volVectorField::New
            (
                IOobject::scopedName(cloudName, "U"),
                IOobject::NO_REGISTER,
                mesh,
                dimensionedVector(dimVelocity, Zero),
                fvPatchFieldBase::zeroGradientType()
            );
            auto& U = tU.ref();

            U.primitiveFieldRef() = uAverage.primitiveField();
            U.correctBoundaryConditions();

            surfaceScalarField phi
            (
                IOobject::scopedName(cloudName, "phi"),
                linearInterpolate(U) & mesh.Sf()
            );

            if (applyGravity_)
            {
                phiCorrect_.ref() -= phiGByA();
            }

            forAll(phiCorrect_(), facei)
            {
                // Current and correction fluxes
                const scalar phiCurr = phi[facei];
                scalar& phiCorr = phiCorrect_.ref()[facei];

                // Don't limit if the correction is in the opposite direction to
                // the flux. We need all the help we can get in this state.
                if (phiCurr*phiCorr < 0)
                {}

                // If the correction and the flux are in the same direction then
                // don't apply any more correction than is already present in
                // the flux.
                else if (phiCorr > 0)
                {
                    phiCorr = max(phiCorr - phiCurr, 0);
                }
                else
                {
                    phiCorr = min(phiCorr - phiCurr, 0);
                }
            }

            if (applyGravity_)
            {
                phiCorrect_.ref() += phiGByA();
            }
        }

        // correction velocity
        uCorrect_ = volVectorField::New
        (
            IOobject::scopedName(cloudName, "uCorrect"),
            IOobject::NO_REGISTER,
            fvc::reconstruct(phiCorrect_())
        );
        uCorrect_->correctBoundaryConditions();

        //Info << endl;
        //Info << "     alpha: " << alpha_.primitiveField() << endl;
        //Info << "phiCorrect: " << phiCorrect_->primitiveField() << endl;
        //Info << "  uCorrect: " << uCorrect_->primitiveField() << endl;
        //Info << endl;
    }
    else
    {
        alpha_.oldTime();
        phiCorrect_.clear();
        uCorrect_.clear();
    }
}


template<class CloudType>
Foam::vector Foam::PackingModels::Implicit<CloudType>::velocityCorrection
(
    typename CloudType::parcelType& p,
    const scalar deltaT
) const
{
    const fvMesh& mesh = this->owner().mesh();

    // containing tetrahedron and parcel coordinates within
    const label celli = p.cell();
    const label facei = p.tetFace();

    // cell velocity
    const vector U = uCorrect_()[celli];

    // face geometry
    vector nHat = mesh.faces()[facei].areaNormal(mesh.points());
    const scalar nMag = mag(nHat);
    nHat /= nMag;

    // get face flux
    scalar phi;
    const label patchi = mesh.boundaryMesh().whichPatch(facei);
    if (patchi == -1)
    {
        phi = phiCorrect_()[facei];
    }
    else
    {
        phi =
            phiCorrect_().boundaryField()[patchi]
            [
                mesh.boundaryMesh()[patchi].whichFace(facei)
            ];
    }

    // interpolant equal to 1 at the cell centre and 0 at the face
    const scalar t = p.coordinates()[0];

    // the normal component of the velocity correction is interpolated linearly
    // the tangential component is equal to that at the cell centre
    return U + (1.0 - t)*nHat*(phi/nMag - (U & nHat));
}


// ************************************************************************* //
