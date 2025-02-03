/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020,2025 OpenCFD Ltd.
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

#include "polyMesh.H"
#include "transform.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Update this with w2 if w2 nearer to pt.
template<class TrackingData>
inline bool Foam::wallPoint::update
(
    const wallPoint& w2,
    const scalar dist2,
    const scalar tol,
    TrackingData& td
)
{
    if (!valid(td))
    {
        // current not yet set so use any value
        distSqr_ = dist2;
        origin_ = w2.origin();

        return true;
    }

    const scalar diff = distSqr_ - dist2;

    if (diff < 0)
    {
        // already nearer to pt
        return false;
    }

    if ((diff < SMALL) || ((distSqr_ > SMALL) && (diff/distSqr_ < tol)))
    {
        // don't propagate small changes
        return false;
    }
    else
    {
        // update with new values
        distSqr_ = dist2;
        origin_ = w2.origin();

        return true;
    }
}


// Update this with w2 if w2 nearer to pt.
template<class TrackingData>
inline bool Foam::wallPoint::update
(
    const point& pt,
    const wallPoint& w2,
    const scalar tol,
    TrackingData& td
)
{
    return update
    (
        w2,
        magSqr(pt - w2.origin()),
        tol,
        td
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::wallPoint::wallPoint()
:
    origin_(point::max),
    distSqr_(-GREAT)
{}


inline Foam::wallPoint::wallPoint(const point& origin, const scalar distSqr)
:
    origin_(origin),
    distSqr_(distSqr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TrackingData>
inline bool Foam::wallPoint::valid(TrackingData& td) const
{
    return distSqr_ > -SMALL;
}


// Checks for cyclic faces
template<class TrackingData>
inline bool Foam::wallPoint::sameGeometry
(
    const polyMesh&,
    const wallPoint& w2,
    const scalar tol,
    TrackingData& td
) const
{
    const scalar diff = mag(distSqr() - w2.distSqr());

    if (diff < SMALL)
    {
        return true;
    }
    else
    {
        if ((distSqr() > SMALL) && ((diff/distSqr()) < tol))
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}


template<class TrackingData>
inline void Foam::wallPoint::leaveDomain
(
    const polyMesh&,
    const polyPatch&,
    const label,
    const point& faceCentre,
    TrackingData& td
)
{
    origin_ -= faceCentre;
}


template<class TrackingData>
inline void Foam::wallPoint::transform
(
    const polyMesh&,
    const tensor& rotTensor,
    TrackingData& td
)
{
    origin_ = Foam::transform(rotTensor, origin_);
}


// Update absolute geometric quantities. Note that distance (distSqr_)
// is not affected by leaving/entering domain.
template<class TrackingData>
inline void Foam::wallPoint::enterDomain
(
    const polyMesh&,
    const polyPatch&,
    const label,
    const point& faceCentre,
    TrackingData& td
)
{
    // back to absolute form
    origin_ += faceCentre;
}


// Update this with w2 if w2 nearer to pt.
template<class TrackingData>
inline bool Foam::wallPoint::updateCell
(
    const polyMesh& mesh,
    const label thisCelli,
    const label neighbourFacei,
    const wallPoint& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    return
        update
        (
            mesh.cellCentres()[thisCelli],
            neighbourWallInfo,
            tol,
            td
        );
}


// Update this with w2 if w2 nearer to pt.
template<class TrackingData>
inline bool Foam::wallPoint::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const label neighbourCelli,
    const wallPoint& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    return
        update
        (
            mesh.faceCentres()[thisFacei],
            neighbourWallInfo,
            tol,
            td
        );
}


// Update this with w2 if w2 nearer to pt.
template<class TrackingData>
inline bool Foam::wallPoint::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const wallPoint& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    return
        update
        (
            mesh.faceCentres()[thisFacei],
            neighbourWallInfo,
            tol,
            td
        );
}


template<class TrackingData>
inline bool Foam::wallPoint::equal
(
    const wallPoint& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


template<class TrackingData>
inline bool Foam::wallPoint::interpolate
(
    const polyMesh&,
    const point& pt,
    const label i0,
    const wallPoint& f0,
    const label i1,
    const wallPoint& f1,
    const scalar weight,
    const scalar tol,
    TrackingData& td
)
{
    if (f0.valid(td))
    {
        if (f1.valid(td))
        {
            // What do we do about the origin? Should be consistent with the
            // distance? But then interpolating between a point 'on the left'
            // and a point 'on the right' the interpolate might just be on
            // top of us so suddenly setting the distance to be zero...

            const scalar dist2 =
                sqr(lerp(sqrt(f0.distSqr()), sqrt(f1.distSqr()), weight));

            if (f0.distSqr() < f1.distSqr())
            {
                // Update with interpolated distance and f0 origin
                return update(f0, dist2, tol, td);
            }
            else
            {
                return update(f1, dist2, tol, td);
            }
        }
        else
        {
            return update(f0, f0.distSqr(), tol, td);
        }
    }
    else if (f1.valid(td))
    {
        return update(f1, f1.distSqr(), tol, td);
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::wallPoint::operator==
(
    const wallPoint& rhs
) const
{
    return origin_ == rhs.origin_;
}


inline bool Foam::wallPoint::operator!=
(
    const wallPoint& rhs
) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
