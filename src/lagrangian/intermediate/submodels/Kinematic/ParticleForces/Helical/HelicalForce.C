/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "HelicalForce.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::HelicalForce<CloudType>::HelicalForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    axisPointPtr_(Function1<vector>::New("axisPoint", this->coeffs(), &mesh)),
    axisDirPtr_(Function1<vector>::New("axisDir", this->coeffs(), &mesh)),
    omegaPtr_(Function1<scalar>::New("omega", this->coeffs(), &mesh)),
    UaxialPtr_
    (
        Function1<scalar>::NewIfPresent("Uaxial", this->coeffs(), &mesh)
    ),
    tauRotationalPtr_
    (
        Function1<scalar>::NewIfPresent("tauRotational", this->coeffs(), &mesh)
    ),
    tauAxialPtr_
    (
        Function1<scalar>::NewIfPresent("tauAxial", this->coeffs(), &mesh)
    ),
    radialDampingPtr_
    (
        Function1<scalar>::NewIfPresent("radialDamping", this->coeffs(), &mesh)
    )
{
    if (tauAxialPtr_ && !UaxialPtr_)
    {
        FatalIOErrorInFunction(this->coeffs())
            << "If tauAxial is specified, Uaxial must also be specified." << nl
            << error(FatalIOError);
    }
}


template<class CloudType>
Foam::HelicalForce<CloudType>::HelicalForce
(
    const HelicalForce<CloudType>& df
)
:
    ParticleForce<CloudType>(df),
    axisPointPtr_(df.axisPointPtr_.clone()),
    axisDirPtr_(df.axisDirPtr_.clone()),
    omegaPtr_(df.omegaPtr_.clone()),
    UaxialPtr_(df.UaxialPtr_.clone()),
    tauRotationalPtr_(df.tauRotationalPtr_.clone()),
    tauAxialPtr_(df.tauAxialPtr_.clone()),
    radialDampingPtr_(df.radialDampingPtr_.clone())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::HelicalForce<CloudType>::~HelicalForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::HelicalForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    const scalar t = this->mesh().time().timeOutputValue();

    const vector e(normalised(axisDirPtr_->value(t)));

    if (mag(e) < VSMALL || !std::isfinite(mag(e)))
    {
        FatalIOErrorInFunction(this->coeffs())
            << "The axisDir vector cannot be zero" << nl
            << "axisDir: " << axisDirPtr_->value(t) << nl
            << "time: " << t << nl
            << error(FatalIOError);
    }

    const vector& X = p.position();
    const vector& U = p.U();


    // Calculate the radial distance between the parcel and rotational axis [m]
    const vector axisPoint(axisPointPtr_->value(t));
    const scalar s = e & (X - axisPoint);  // axial projection, s=e . (x-r0)
    const vector Pax = axisPoint + s*e;    // closest point on axis P(x)
    vector R = X - Pax;                    // radial vector
    scalar r = mag(R);                     // radius, r =|R|


    // Calculate the unit vector of the radial vector, Rhat
    vector Rhat(vector::zero);
    if (r > VSMALL)  // if parcel is farther than VSMALL from the rotation axis
    {
        Rhat = R/r;
    }
    else  // if parcel is effectively on the rotation axis, pick something
    {
        const vector n
        (
            (mag(e.x()) < 0.9) ? vector(1,0,0) : vector(0,1,0)
        );
        // When the parcel is on the axis, use the cross product to generate a
        // perpendicular direction, ensuring a valid radial vector for further
        // calculations
        R = e^n;
        r = mag(R);
        Rhat = R/r;
    }


    // Calculate the velocity components of the parcel
    const scalar Uaxial = (U & e);
    const vector Utransverse(U - Uaxial*e);
    const scalar Uradial = (Utransverse & Rhat);
    const vector Utangent(Utransverse - Uradial*Rhat);


    // Accumulate the acceleration components for the parcel
    const scalar omega = omegaPtr_->value(t);
    // Centripetal acceleration: directed inward toward the axis of rotation
    // The negative sign ensures the acceleration points from the particle's
    // position toward the axis.
    vector a(-sqr(omega)*R);


    if (tauRotationalPtr_)  // try to adjust for the desired tangential velocity
    {
        const scalar tau = tauRotationalPtr_->value(t);

        if (tau < SMALL || !std::isfinite(tau))
        {
            FatalIOErrorInFunction(this->coeffs())
                << "tauRotational cannot be zero or negative" << nl
                << "tauRotational: " << tau << nl
                << "time: " << t << nl
                << error(FatalIOError);
        }

        const vector UtangentDesired(omega*(e^R));
        a += (UtangentDesired - Utangent)/tau;
    }


    if (tauAxialPtr_ )  // try to adjust for desired axial velocity
    {
        if (!UaxialPtr_)
        {
            FatalIOErrorInFunction(this->coeffs())
                << "If tauAxial is specified, Uaxial must also be specified."
                << exit(FatalIOError);
        }

        const scalar tau = tauAxialPtr_->value(t);

        if (tau < SMALL || !std::isfinite(tau))
        {
            FatalIOErrorInFunction(this->coeffs())
                << "tauAxial cannot be zero or negative" << nl
                << "tauAxial: " << tau << nl
                << "time: " << t << nl
                << error(FatalIOError);
        }

        const scalar UaxialDesired = UaxialPtr_->value(t);
        a += ((UaxialDesired - Uaxial)/tau)*e;
    }


    // Apply radial damping to suppress radial oscillations
    if (radialDampingPtr_)
    {
        const scalar damper = radialDampingPtr_->value(t);

        if (damper < SMALL || !std::isfinite(damper))
        {
            FatalIOErrorInFunction(this->coeffs())
                << "radialDamping cannot be zero or negative" << nl
                << "radialDamping: " << damper << nl
                << "time: " << t << nl
                << error(FatalIOError);
        }

        a += (-damper*Uradial)*Rhat;
    }


    return forceSuSp(mass*a, Zero);  // (Sp, Su)
}


// ************************************************************************* //
