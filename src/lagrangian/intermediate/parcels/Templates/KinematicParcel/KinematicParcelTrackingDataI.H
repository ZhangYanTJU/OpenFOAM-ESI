/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "AveragingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
inline Foam::KinematicParcel<ParcelType>::trackingData::trackingData
(
    const TrackCloudType& cloud,
    trackPart part
)
:
    ParcelType::trackingData(cloud),
    rhoInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.rho()
        )
    ),
    UInterp_
    (
        interpolation<vector>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.U()
        )
    ),
    muInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.mu()
        )
    ),
    rhoc_(Zero),
    Uc_(Zero),
    muc_(Zero),

    volumeAverage_
    (
        AveragingMethod<scalar>::New
        (
            IOobject
            (
                IOobject::scopedName(cloud.name(), "volumeAverage"),
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.solution().dict(),
            cloud.mesh()
        )
    ),
    radiusAverage_
    (
        AveragingMethod<scalar>::New
        (
            IOobject
            (
                IOobject::scopedName(cloud.name(), "radiusAverage"),
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.solution().dict(),
            cloud.mesh()
        )
    ),
    rhoAverage_
    (
        AveragingMethod<scalar>::New
        (
            IOobject
            (
                IOobject::scopedName(cloud.name(), "rhoAverage"),
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.solution().dict(),
            cloud.mesh()
        )
    ),
    uAverage_
    (
        AveragingMethod<vector>::New
        (
            IOobject
            (
                IOobject::scopedName(cloud.name(), "uAverage"),
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.solution().dict(),
            cloud.mesh()
        )
    ),
    uSqrAverage_
    (
        AveragingMethod<scalar>::New
        (
            IOobject
            (
                IOobject::scopedName(cloud.name(), "uSqrAverage"),
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.solution().dict(),
            cloud.mesh()
        )
    ),
    frequencyAverage_
    (
        AveragingMethod<scalar>::New
        (
            IOobject
            (
                IOobject::scopedName(cloud.name(), "frequencyAverage"),
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.solution().dict(),
            cloud.mesh()
        )
    ),
    massAverage_
    (
        AveragingMethod<scalar>::New
        (
            IOobject
            (
                IOobject::scopedName(cloud.name(), "massAverage"),
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.solution().dict(),
            cloud.mesh()
        )
    ),

    g_(cloud.g().value()),
    part_(part)
{}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::KinematicParcel<ParcelType>::trackingData::rhoInterp() const
{
    return *rhoInterp_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::vector>&
Foam::KinematicParcel<ParcelType>::trackingData::UInterp() const
{
    return *UInterp_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::KinematicParcel<ParcelType>::trackingData::muInterp() const
{
    return *muInterp_;
}


template<class ParcelType>
inline const Foam::vector&
Foam::KinematicParcel<ParcelType>::trackingData::g() const
{
    return g_;
}


template<class ParcelType>
inline Foam::scalar
Foam::KinematicParcel<ParcelType>::trackingData::rhoc() const
{
    return rhoc_;
}


template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::trackingData::rhoc()
{
    return rhoc_;
}


template<class ParcelType>
inline const Foam::vector&
Foam::KinematicParcel<ParcelType>::trackingData::Uc() const
{
    return Uc_;
}


template<class ParcelType>
inline Foam::vector& Foam::KinematicParcel<ParcelType>::trackingData::Uc()
{
    return Uc_;
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicParcel<ParcelType>::trackingData::muc() const
{
    return muc_;
}


template<class ParcelType>
inline Foam::scalar& Foam::KinematicParcel<ParcelType>::trackingData::muc()
{
    return muc_;
}


template<class ParcelType>
inline typename Foam::KinematicParcel<ParcelType>::trackingData::trackPart
Foam::KinematicParcel<ParcelType>::trackingData::part() const
{
    return part_;
}


template<class ParcelType>
inline typename Foam::KinematicParcel<ParcelType>::trackingData::trackPart&
Foam::KinematicParcel<ParcelType>::trackingData::part()
{
    return part_;
}


template<class ParcelType>
template<class TrackCloudType>
inline void Foam::KinematicParcel<ParcelType>::trackingData::
updateAverages
(
    const TrackCloudType& cloud
)
{
    // zero the sums
    volumeAverage_() = 0;
    radiusAverage_() = 0;
    rhoAverage_() = 0;
    uAverage_() = Zero;
    uSqrAverage_() = 0;
    frequencyAverage_() = 0;
    massAverage_() = 0;

    // temporary weights
    autoPtr<AveragingMethod<scalar>> weightAveragePtr
    (
        AveragingMethod<scalar>::New
        (
            IOobject
            (
                IOobject::scopedName(cloud.name(), "weightAverage"),
                cloud.db().time().timeName(),
                cloud.mesh()
            ),
            cloud.solution().dict(),
            cloud.mesh()
        )
    );
    AveragingMethod<scalar>& weightAverage = weightAveragePtr();

    // averaging sums
    for (const typename TrackCloudType::parcelType& p : cloud)
    {
        const tetIndices tetIs = p.currentTetIndices();

        const scalar m = p.nParticle()*p.mass();

        volumeAverage_->add(p.coordinates(), tetIs, p.nParticle()*p.volume());
        rhoAverage_->add(p.coordinates(), tetIs, m*p.rho());
        uAverage_->add(p.coordinates(), tetIs, m*p.U());
        massAverage_->add(p.coordinates(), tetIs, m);
    }
    volumeAverage_->average();
    massAverage_->average();
    rhoAverage_->average(*massAverage_);
    uAverage_->average(*massAverage_);

    // squared velocity deviation
    for (const typename TrackCloudType::parcelType& p : cloud)
    {
        const tetIndices tetIs = p.currentTetIndices();

        const vector u = uAverage_->interpolate(p.coordinates(), tetIs);

        uSqrAverage_->add
        (
            p.coordinates(),
            tetIs,
            p.nParticle()*p.mass()*magSqr(p.U() - u)
        );
    }
    uSqrAverage_->average(*massAverage_);

    // sauter mean radius
    radiusAverage_() = volumeAverage_();
    weightAverage = 0;
    for (const typename TrackCloudType::parcelType& p : cloud)
    {
        const tetIndices tetIs = p.currentTetIndices();

        weightAverage.add
        (
            p.coordinates(),
            tetIs,
            p.nParticle()*pow(p.volume(), 2.0/3.0)
        );
    }
    weightAverage.average();
    radiusAverage_->average(weightAverage);

    // collision frequency
    weightAverage = 0;
    for (const typename TrackCloudType::parcelType& p : cloud)
    {
        const tetIndices tetIs = p.currentTetIndices();

        const scalar a = volumeAverage_->interpolate(p.coordinates(), tetIs);
        const scalar r = radiusAverage_->interpolate(p.coordinates(), tetIs);
        const vector u = uAverage_->interpolate(p.coordinates(), tetIs);

        const scalar f = 0.75*a/pow3(r)*sqr(0.5*p.d() + r)*mag(p.U() - u);

        frequencyAverage_->add(p.coordinates(), tetIs, p.nParticle()*f*f);

        weightAverage.add(p.coordinates(), tetIs, p.nParticle()*f);
    }
    frequencyAverage_->average(weightAverage);
}


// ************************************************************************* //
