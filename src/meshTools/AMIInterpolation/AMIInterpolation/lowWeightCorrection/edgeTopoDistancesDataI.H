/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "polyMesh.H"
#include "transform.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class PrimitivePatchType>
inline
Foam::edgeTopoDistancesData<Type, PrimitivePatchType>::edgeTopoDistancesData()
:
    distance_(),
    data_()
{}


template<class Type, class PrimitivePatchType>
inline
Foam::edgeTopoDistancesData<Type, PrimitivePatchType>::edgeTopoDistancesData
(
    const labelList& distance,
    const List<Type>& data
)
:
    distance_(distance),
    data_(data)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class PrimitivePatchType>
template<class TrackingData>
inline bool Foam::edgeTopoDistancesData<Type, PrimitivePatchType>::valid
(
    TrackingData& td
) const
{
    return distance_.size();
}


template<class Type, class PrimitivePatchType>
template<class TrackingData>
inline void Foam::edgeTopoDistancesData<Type, PrimitivePatchType>::transform
(
    const polyMesh& mesh,
    const PrimitivePatchType& patch,
    const tensor& rotTensor,
    const scalar tol,
    TrackingData& td
)
{}


template<class Type, class PrimitivePatchType>
template<class TrackingData>
inline bool Foam::edgeTopoDistancesData<Type, PrimitivePatchType>::updateEdge
(
    const polyMesh& mesh,
    const PrimitivePatchType& patch,
    const label edgeI,
    const label facei,
    const edgeTopoDistancesData<Type, PrimitivePatchType>& faceInfo,
    const scalar tol,
    TrackingData& td
)
{
    // From face to edge

    if (distance_.size() >= td.n_)
    {
        return false;
    }

    const label oldSize = data_.size();

    forAll(faceInfo.data_, i)
    {
        const auto& d = faceInfo.data_[i];

        if (!data_.found(d))
        {
            data_.append(d);
            distance_.append(faceInfo.distance_[i] + 1);
        }
    }
    return data_.size() > oldSize;
}


template<class Type, class PrimitivePatchType>
template<class TrackingData>
inline bool Foam::edgeTopoDistancesData<Type, PrimitivePatchType>::updateEdge
(
    const polyMesh& mesh,
    const PrimitivePatchType& patch,
    const edgeTopoDistancesData<Type, PrimitivePatchType>& edgeInfo,
    const bool sameOrientation,
    const scalar tol,
    TrackingData& td
)
{
    // From edge to edge (e.g. coupled edges)

    if (distance_.size() >= td.n_)
    {
        return false;
    }

    const label oldSize = data_.size();

    forAll(edgeInfo.data_, i)
    {
        const auto& d = edgeInfo.data_[i];

        if (!data_.found(d))
        {
            data_.append(d);
            distance_.append(edgeInfo.distance_[i]);
        }
    }

    return data_.size() > oldSize;
}


template<class Type, class PrimitivePatchType>
template<class TrackingData>
inline bool Foam::edgeTopoDistancesData<Type, PrimitivePatchType>::updateFace
(
    const polyMesh& mesh,
    const PrimitivePatchType& patch,
    const label facei,
    const label edgeI,
    const edgeTopoDistancesData<Type, PrimitivePatchType>& edgeInfo,
    const scalar tol,
    TrackingData& td
)
{
    // From edge to edge (e.g. coupled edges)

    if (distance_.size() >= td.n_)
    {
        return false;
    }

    const label oldSize = data_.size();

    forAll(edgeInfo.data_, i)
    {
        const auto& d = edgeInfo.data_[i];

        if (!data_.found(d))
        {
            data_.append(d);
            distance_.append(edgeInfo.distance_[i]);
        }
    }

    return data_.size() > oldSize;
}


template<class Type, class PrimitivePatchType>
template<class TrackingData>
inline bool Foam::edgeTopoDistancesData<Type, PrimitivePatchType>::equal
(
    const edgeTopoDistancesData<Type, PrimitivePatchType>& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type, class PrimitivePatchType>
inline bool Foam::edgeTopoDistancesData<Type, PrimitivePatchType>::operator==
(
    const Foam::edgeTopoDistancesData<Type, PrimitivePatchType>& rhs
) const
{
    return distance() == rhs.distance() && data() == rhs.data();
}


template<class Type, class PrimitivePatchType>
inline bool Foam::edgeTopoDistancesData<Type, PrimitivePatchType>::operator!=
(
    const Foam::edgeTopoDistancesData<Type, PrimitivePatchType>& rhs
) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
