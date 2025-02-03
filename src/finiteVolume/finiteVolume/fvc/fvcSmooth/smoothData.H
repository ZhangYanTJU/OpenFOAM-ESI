/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

Class
    Foam::smoothData

Description
    Helper class used by the fvc::smooth and fvc::spread functions.

SourceFiles
    smoothData.H
    smoothDataI.H

\*---------------------------------------------------------------------------*/

#ifndef smoothData_H
#define smoothData_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class smoothData;

Ostream& operator<<(Ostream&, const smoothData&);
Istream& operator>>(Istream&, smoothData&);

/*---------------------------------------------------------------------------*\
                       Class smoothData Declaration
\*---------------------------------------------------------------------------*/

class smoothData
{
    // Private Data

        scalar value_;


    // Private Member Functions

        //- Update gets information from neighbouring face/cell and
        //- uses this to update itself (if necessary) and return true
        template<class TrackingData>
        inline bool update
        (
            const smoothData& svf,
            const scalar scale,
            const scalar tol,
            TrackingData& td
        );


public:

    //- Class used to pass additional data in
    class trackData
    {
    public:

        //- Cut off distance
        scalar maxRatio;
    };


    // Generated Methods

        //- Copy construct
        smoothData(const smoothData&) = default;

        //- Copy assignment
        smoothData& operator=(const smoothData&) = default;


    // Constructors

        //- Default construct
        inline smoothData();

        //- Construct from inverse field value
        inline smoothData(const scalar value);


    // Member Functions

        // Access

            //- Return value
            scalar value() const
            {
                return value_;
            }

            //- The value
            scalar& value()
            {
                return value_;
            }


        // Needed by FaceCellWave

            //- Changed or contains original (invalid) value
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Check for identical geometrical data (eg, cyclics checking)
            template<class TrackingData>
            inline bool sameGeometry
            (
                const polyMesh&,
                const smoothData&,
                const scalar,
                TrackingData& td
            ) const;

            //- Convert any absolute coordinates into relative to
            //- (patch)face centre
            template<class TrackingData>
            inline void leaveDomain
            (
                const polyMesh&,
                const polyPatch&,
                const label patchFacei,
                const point& faceCentre,
                TrackingData& td
            );

            //- Reverse of leaveDomain
            template<class TrackingData>
            inline void enterDomain
            (
                const polyMesh&,
                const polyPatch&,
                const label patchFacei,
                const point& faceCentre,
                TrackingData& td
            );

            //- Apply rotation matrix to any coordinates
            template<class TrackingData>
            inline void transform
            (
                const polyMesh&,
                const tensor&,
                TrackingData& td
            );

            //- Influence of neighbouring face
            template<class TrackingData>
            inline bool updateCell
            (
                const polyMesh&,
                const label thisCelli,
                const label neighbourFacei,
                const smoothData& svf,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of neighbouring cell
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const label neighbourCelli,
                const smoothData& svf,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const smoothData& svf,
                const scalar tol,
                TrackingData& td
            );

            //- Test for equality, with TrackingData
            template<class TrackingData>
            inline bool equal(const smoothData&, TrackingData& td) const;

            //- Interpolate between two values (lerp). Returns true if
            //- causes changes. Not sure if needs to be specialised between
            //- face and cell and what index is needed...
            template<class TrackingData>
            inline bool interpolate
            (
                const polyMesh&,
                const point& pt,
                const label i0,
                const smoothData& f0,
                const label i1,
                const smoothData& f1,
                const scalar weight,
                const scalar tol,
                TrackingData& td
            );


    // Member Operators

        //- Assign new value
        void operator=(const scalar value)
        {
            value_ = value;
        }

        //- Test for equality
        inline bool operator==(const smoothData&) const;

        //- Test for inequality
        inline bool operator!=(const smoothData&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream& os, const smoothData& rhs)
        {
            return os << rhs.value_;
        }

        friend Istream& operator>>(Istream& is, smoothData& rhs)
        {
            return is >> rhs.value_;
        }
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for smoothData
template<> struct is_contiguous<smoothData> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "smoothDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
