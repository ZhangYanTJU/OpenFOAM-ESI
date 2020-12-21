/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Class
    Foam::PDRfitMeshScan

Description
    Scanning of obstacles in a single direction for PDRfitMesh

SourceFiles
    PDRfitMeshScan.C

\*---------------------------------------------------------------------------*/

#ifndef PDRfitMeshScan_H
#define PDRfitMeshScan_H

#include "scalarList.H"
#include "DynamicList.H"
#include "MinMax.H"
#include "PDRblock.H"
#include "PDRfitMeshParams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class PDRfitMeshScan Declaration
\*---------------------------------------------------------------------------*/

class PDRfitMeshScan
{
    // Private Data

        //- Lower and upper limits
        scalarMinMax limits_;

        //- Absolute lower limit (eg, for the ground plane)
        scalar hard_min_ = -GREAT;

        //- Number of steps
        label nsteps_ = 0;

        //- Step width = limits / nsteps
        scalar stepWidth_ = 0;

        //- Minimum per cell face area. Defined via PDRfitMeshParams
        scalar minFaceArea_ = 0;

        //- Area-weighted position
        DynamicList<scalar> weightedPos_;

        //- Area-weights for position
        DynamicList<scalar> totalArea_;

        //- Area-weighted half-cell position
        DynamicList<scalar> weightedPos2_;

        //- Area-weights for half-cell position
        DynamicList<scalar> totalArea2_;


    // Private Member Functions

        //- Area-averaged position at index
        inline scalar averagePosition(const label i) const
        {
            return
            (
                totalArea_[i] < VSMALL
              ? 0
              : (weightedPos_[i] / totalArea_[i])
            );
        }

        //- Area-averaged half-position at index
        inline scalar averageHalfPosition(const label i) const
        {
            return
            (
                totalArea2_[i] < VSMALL
              ? 0
              : (weightedPos2_[i] / totalArea2_[i])
            );
        }


        //- Add coordinate and weight. Already sanity checked
        void addCoord(const scalar coord, const scalar weight);


public:

    // Public Data

        //- Output verbosity
        static bool verbose_;

        //- Use gridControl as per PDRblock, with EXPAND_RELATIVE
        typedef PDRblock::gridControl gridControl;


    // Constructors

        //- Default construct
        PDRfitMeshScan()
        :
            limits_(GREAT, -GREAT),
            hard_min_(-GREAT),
            nsteps_(0),
            stepWidth_(0)
        {}


    // Static Member Functions

        //- Get verbosity
        static inline bool verbose() noexcept
        {
            return verbose_;
        }

        //- Set verbosity
        static inline bool verbose(const bool val) noexcept
        {
            bool old(verbose_);
            verbose_ = val;
            return old;
        }

        //- The expansion name
        static word expansionName();


    // Member Functions

        //- The min/max limits
        const scalarMinMax& limits() const noexcept
        {
            return limits_;
        }

        //- The min/max limits
        scalarMinMax& limits() noexcept
        {
            return limits_;
        }

        //- The step-width
        scalar stepWidth() const noexcept
        {
            return stepWidth_;
        }

        //- Set the hard-min
        void hard_min(const scalar val) noexcept
        {
            hard_min_ = val;
        }


        //- Reset to initial state
        void reset();

        //- Adjust limits to accommodate the specified point coordinates
        void adjustLimits(const scalar point0, const scalar point1);

        //- Define number of steps (and step-size) according to the
        //- current limits and the specified parameters
        void resize(const scalar cellWidth, const PDRfitMeshParams& pars);

        void addArea
        (
            const scalar weight,
            const scalar point0,
            const scalar point1
        );

        //- Calculate equivalent grid control from the scanned planes
        PDRblock::gridControl
        calcGridControl
        (
            const scalar cellWidth,
            const scalar max_zone,
            const PDRfitMeshParams& pars
        );

        void print(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
