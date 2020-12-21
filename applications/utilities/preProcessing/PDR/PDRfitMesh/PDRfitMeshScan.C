/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 Shell Research Ltd.
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

#include "PDRfitMeshScan.H"
#include "PDRfitMeshParams.H"
#include "PDRblock.H"
#include "PDRobstacle.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::PDRfitMeshScan::verbose_ = false;


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::word Foam::PDRfitMeshScan::expansionName()
{
    return PDRblock::expansionNames_[PDRblock::EXPAND_RELATIVE];
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::PDRfitMeshScan::addCoord
(
    const scalar coord,
    const scalar weight
)
{
    if (!limits_.contains(coord))
    {
        return;
    }

    const scalar coordOffset = (coord - limits_.min());

    if (coordOffset < 0)
    {
        FatalErrorInFunction
            << "Negative coordinate! " << coordOffset << nl
            << exit(FatalError);
    }

    // We divide the range into a series of steps, and decide
    // which step this coord is in


    // Note: would be simpler to record the totals only in the
    // half steps in this routine then get the full-step values
    // at the end by adding half steps


    // Add the area to the total already found in this step
    // and update the weighted average position of these obstacle faces

    const label stepi1 = 2*floor(coordOffset / stepWidth_);


    // Also do half-step offset in case several faces
    // are around the step boundary.
    //
    // Stored in odd-numbered elements

    const label stepi2 =
        Foam::max
        (
            0,
            1 + 2*floor((coordOffset - 0.5*stepWidth_) / stepWidth_)
        );


    // Last, keep track of totals in half-step ranges

    const label stepih = floor(coordOffset / stepWidth_);


    const label maxSize = Foam::max(stepi1, stepi2);

    if (weightedPos_.size() <= maxSize)
    {
        weightedPos_.resize(maxSize+1, Zero);
        totalArea_.resize(maxSize+1, Zero);

        weightedPos2_.resize(maxSize+1, Zero);
        totalArea2_.resize(maxSize+1, Zero);
    }

    weightedPos_[stepi1] += (coord * weight);
    totalArea_[stepi1] += weight;

    weightedPos_[stepi2] += (coord * weight);
    totalArea_[stepi2] += weight;

    weightedPos2_[stepih] += (coord * weight);
    totalArea2_[stepih] += weight;
}


// void Foam::PDRfitMeshScan::clear()
// {
//     weightedPos_.clear();  totalArea_.clear();
//     weightedPos2_.clear(); totalArea2_.clear();
// }


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PDRfitMeshScan::reset()
{
    limits_.min() = GREAT;
    limits_.max() = -GREAT;

    hard_min_ = -GREAT;

    nsteps_ = 0;
    stepWidth_ = 0;
    minFaceArea_ = 0;

    weightedPos_.clear();
    totalArea_.clear();

    weightedPos2_.clear();
    totalArea2_.clear();
}


void Foam::PDRfitMeshScan::resize
(
    const scalar cellWidth,
    const PDRfitMeshParams& pars
)
{
    if (limits_.valid())
    {
        nsteps_ =
        (
            limits_.mag() / (mag(cellWidth) * pars.areaWidthFactor)
          + 0.5
        );

        nsteps_ = max(nsteps_, pars.nCellsMin);

        stepWidth_ = limits_.mag() / nsteps_;

        minFaceArea_ = pars.minFaceArea;
    }
    else
    {
        nsteps_ = 0;
        stepWidth_ = 0;
        minFaceArea_ = 0;

        WarningInFunction
            << "No valid limits defined" << endl;
    }

    weightedPos_.clear();
    totalArea_.clear();

    weightedPos2_.clear();
    totalArea2_.clear();

    // resize
    weightedPos_.resize(2*nsteps_+1, Zero);
    totalArea_.resize(2*nsteps_+1, Zero);

    weightedPos2_.resize(2*nsteps_+1, Zero);
    totalArea2_.resize(2*nsteps_+1, Zero);
}


void Foam::PDRfitMeshScan::adjustLimits
(
    const scalar point0,
    const scalar point1
)
{
    limits_.add(point0);
    limits_.add(point1);
}


void Foam::PDRfitMeshScan::addArea
(
    const scalar area,
    const scalar point0,
    const scalar point1
)
{
    if (!nsteps_)
    {
        FatalErrorInFunction
            << "No step-size defined" << nl
            << exit(FatalError);
    }

    if (area < minFaceArea_ * sqr(stepWidth_))
    {
        return;
    }

    addCoord(point0, area);
    addCoord(point1, area);
}


void Foam::PDRfitMeshScan::print(Ostream& os) const
{
    if (nsteps_)
    {
        os  << "steps:" << nsteps_;
    }

    if (limits_.valid())
    {
        os  << " limits:" << limits_;
    }

    scalar totalArea = 0;
    for (const scalar areaValue : totalArea_)
    {
        totalArea += areaValue;
    }

    os  << " area:" << totalArea;
    os  << nl;
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Evaluates (r**n - 1) / (r - 1) including in the limit r=1
inline static scalar rn1r1(const scalar r, const label n)
{
    if (mag(r - 1.0) < SMALL)
    {
        return n;
    }

    return (pow(r, n) - 1.0)/(r - 1.0);
}


inline static scalar end_val
(
    const scalar width,
    const scalar ratio,
    const label n
)
{
    return width / rn1r1(ratio, n);
}


static scalar fit_slope
(
    const scalar left,
    const scalar width,
    const label n,
    scalar& r
)
{
    const scalar wol = width / left;

    // Initial value
    r = (left < width / n) ? 1.01 : 0.99;

    for (int nIter = 0; nIter < 25; ++nIter)
    {
        scalar rm1 = r - 1.0;
        scalar rn = pow(r, n);
        scalar f = (rn - 1.0) * r / rm1 - wol;
        scalar fprime = ((n * rm1 - 1) * rn + 1) / sqr(rm1);

        scalar new_r = r - f / fprime;
        scalar delta = mag(new_r - r);

        r = 0.5 * (r + new_r);

        // InfoErr
        //     << "New r: " << new_r << ' '
        //     << f << ' ' << fprime << ' ' << delta  << ' ' << nIter << nl;

        if (delta <= 1e-3)
        {
            break;
        }
    }

    return end_val(width, 1.0/r, n);
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::PDRblock::gridControl
Foam::PDRfitMeshScan::calcGridControl
(
    const scalar cellWidth,
    const scalar max_zone,
    const PDRfitMeshParams& pars
)
{
    // Prune out small areas
    {
        const scalar minArea = sqr(cellWidth) * pars.minAreaRatio;
        if (verbose())
        {
            Info<< "Checking planes. Min-area: " << minArea << nl;
        }

        for (scalar& areaValue : totalArea_)
        {
            if (areaValue < minArea)
            {
                areaValue = 0;
            }
        }

        for (scalar& areaValue : totalArea2_)
        {
            if (areaValue < minArea)
            {
                areaValue = 0;
            }
        }
    }


    // Initially only the two outer boundaries and absolute limits

    DynamicList<scalar> initialCoord(nsteps_+1);

    initialCoord.resize(4);
    initialCoord[0] = -GREAT;
    initialCoord[1] = Foam::max(limits_.min(), hard_min_);
    initialCoord[2] = limits_.max();
    initialCoord[3] = GREAT;


    // Select largest areas for fitting.
    // - done manually (instead of via sort) to handle full/half-steps
    // - find the largest area, add that to the list and remove its
    //   area (and immediate surroundings) from the input.

    // We have the "full steps overlapping so that we do not miss
    // a group of adjacent faces that are split over a boundary.
    // When we have identified a plane, we need to zero that
    // step, and also the overlapping ones to avoid
    // double-counting. when we have done that twice, we might
    // have completely removed the half-step that was between
    // them. That is why we also have the areas stored in
    // half-steps, so that is not lost.

    while (initialCoord.size()-4 < totalArea_.size())
    {
        label n_big = 0;
        scalar largest = 0;

        // Find largest area, using "real" scanned areas
        // thus checking the stored even locations.

        for (label ii=0; ii < 2*nsteps_; ++ii)
        {
            if (largest < totalArea_[ii])
            {
                largest = totalArea_[ii];
                n_big = ii;
            }
            if (largest < totalArea2_[ii])
            {
                largest = totalArea2_[ii];
                n_big = -ii-1;
            }
        }

        if (mag(largest) < SMALL)
        {
            break;
        }

        if (n_big >= 0)
        {
            // Full step

            const scalar pos = averagePosition(n_big);

            // Remove from future checks
            totalArea_[n_big] = 0;

            if (pos > hard_min_ + 0.4 * cellWidth)
            {
                // Only accept if not close to the ground

                initialCoord.append(pos);

                // Remove this and adjacent areas, which overlap this one
                if (n_big > 1)
                {
                    totalArea_[n_big-1] = 0;
                }
                if (n_big < 2*nsteps_ - 1)
                {
                    totalArea_[n_big+1] = 0;
                }

                // Remove half-steps contained by this step
                totalArea2_[n_big] = 0;
                totalArea2_[n_big+1] = 0;
            }
        }
        else
        {
            // Half step
            n_big = -n_big-1;

            const scalar pos = averageHalfPosition(n_big);

            // Remove from future checks
            totalArea2_[n_big] = 0;

            // On the first passes, this area will be included in full steps,
            // but it may eventually only remain in the half-step,
            // which could be the next largest

            if (pos > hard_min_ + 0.4 * cellWidth)
            {
                // Only accept if not close to the ground

                initialCoord.append(pos);
            }
        }
    }

    // Now sort the positions
    Foam::sort(initialCoord);

    // Avoid small zones near the begin/end positions

    if
    (
        (initialCoord[2] - initialCoord[1])
      < (0.4 * cellWidth)
    )
    {
        // Remove [1] if too close to [2]
        initialCoord.remove(1);
    }

    if
    (
        (
            initialCoord[initialCoord.size()-2]
          - initialCoord[initialCoord.size()-3]
        )
      < (0.4 * cellWidth)
    )
    {
        // Remove [N-3] if too close to [N-2]
        initialCoord.remove(initialCoord.size()-3);
    }


    // This is somewhat questionable (horrible mix of C and C++ addressing).
    //
    // If we have specified a max_zone (max number of cells per segment)
    // we will actually generate more divisions than originally
    // anticipated.

    // So over-dimension arrays by a larger factor (eg, 20-100)
    // which is presumably good enough.

    // FIXME: Needs revisiting (2020-12-16)

    label nCoords = initialCoord.size();

    scalarList coord(50*initialCoord.size(), Zero);
    SubList<scalar>(coord, initialCoord.size()) = initialCoord;

    scalarList width(coord.size(), Zero);
    labelList nDivs(width.size(), Zero);

    //
    // Set up initially uniform zones
    //
    for (label ii=1; ii < (nCoords-2); ++ii)
    {
        width[ii] = coord[ii+1] - coord[ii];
        nDivs[ii] = width[ii] / stepWidth_ * pars.areaWidthFactor + 0.5;

        while (nDivs[ii] < 1)
        {
            // Positions too close - merge them
            --nCoords;

            coord[ii] = 0.5 * (coord[ii+1] + coord[ii]);
            width[ii-1] = coord[ii] - coord[ii-1];
            nDivs[ii-1] =
                width[ii-1] / stepWidth_ * pars.areaWidthFactor + 0.5;

            if (ii < (nCoords - 2))
            {
                for (label iii=ii+1; iii < nCoords-1; ++iii)
                {
                    coord[iii] = coord[iii+1];
                }
                width[ii] = coord[ii+1] - coord[ii];

                nDivs[ii] =
                    width[ii] / stepWidth_ * pars.areaWidthFactor + 0.5;
            }
            else
            {
                nDivs[ii] = 1; // Dummy value to terminate while
            }
        }
    }


    // Initially all steps have equal width (ratio == 1)

    scalarList first(width.size(), Zero);

    for (label ii=1; ii < (nCoords-2); ++ii)
    {
        first[ii] = width[ii] / nDivs[ii];
    }

    scalarList last(first);
    scalarList ratio(width.size(), scalar(1));

    // For diagnostics
    charList marker(width.size(), char('x'));
    marker[0] = '0';


    // Now adjust the inter-cell ratios to reduce cell-size changes
    // between zones

    if (nCoords > 5)
    {
        // The outer zones can be adjusted to fit the adjacent steps

        last[1] = first[2];
        first[nCoords-3] = last[nCoords-4];

        // We adjust the cell size ratio in each zone so that the
        // cell sizes at the ends of the zones fit as well as
        // possible with their neighbours. Work forward through the
        // zones and then repeat a few times so that the adjustment
        // settles down.

        for (int nIter = 0; nIter < 2; ++nIter)
        {
            for (label ii=1; ii < (nCoords-2); ++ii)
            {
                if (nDivs[ii] == 1)
                {
                    // If only one step, no grading is possible

                    marker[ii] = '1';   // For Diagnostics
                    continue;
                }
                marker[ii] = 'L';   // For Diagnostics

                // Try making the in-zone steps equal to the step at
                // the left end to the last of the previous zone.

                // If the step at the right is less than this, then
                // good... (fF this is the penultimate zone, we do
                // not need to worry about the step on the right
                // because the extent of the last zone can be
                // adjusted later to fit.)

                scalar r_fit = 1;

                if
                (
                    (
                        (ii > 1)
                     && mag
                        (
                            log
                            (
                                fit_slope
                                (
                                    last[ii-1],
                                    width[ii],
                                    nDivs[ii],
                                    r_fit
                                )
                            / first[ii+1]
                            )
                        )
                      < mag(log(r_fit))
                    )
                 || (ii == nCoords-3)
                )
                {
                    ratio[ii] = r_fit;
                }
                else
                {
                    marker[ii] = 'R';

                    // otherwise try making the in-zone steps equal to
                    // the step at the right end to the last of the
                    // previous zone. If the step at the right is less
                    // thhan this, then good... (If this is the second
                    // zone, we do not need to worry about the step on
                    // the left because the extent of the first zone
                    // can be adjusted later to fit.)

                    if
                    (
                        mag
                        (
                            log
                            (
                                fit_slope
                                (
                                    first[ii+1],
                                    width[ii],
                                    nDivs[ii],
                                    r_fit
                                )
                              / last[ii-1]
                            )
                        )
                      < mag(log(r_fit))

                     || (ii == 1)
                    )
                    {
                        ratio[ii] = 1.0 / r_fit;
                    }
                    else
                    {
                        ratio[ii] =
                            pow(first[ii+1] / last[ii-1], 1.0/(nDivs[ii]-1));

                        marker[ii] = 'M';
                    }
                }

                first[ii] = end_val(width[ii],     ratio[ii], nDivs[ii]);
                last[ii]  = end_val(width[ii], 1.0/ratio[ii], nDivs[ii]);
            }
        };


        // Adjust the ratio and the width in the first zone to
        // grade from cellWidth at the outside to the first cell of
        // the next zone. Increase the number of steps to keep edge
        // at least two (nEdgeLayers) cell widths from min.
    }


    if (coord[1] > hard_min_ + first[1] / pars.outerRatio)
    {
        nDivs[0] = pars.nEdgeLayers;
        ratio[0] = 1.0 / pars.outerRatio;
        coord[0] = coord[1] - first[1] * rn1r1(pars.outerRatio, nDivs[0]);

        // Reduce the number of steps if we have extended below hard_min_
        while (coord[0] < hard_min_ && nDivs[0] > 1)
        {
            --nDivs[0];
            coord[0] = coord[1] - first[1] * rn1r1(pars.outerRatio, nDivs[0]);
        }

        if (nDivs[0] < label(pars.nEdgeLayers))
        {
            // We had to adjust for hard_min_, so now fit exaxtly
            coord[0] = hard_min_;
        }
    }
    else
    {
        // No room for outer zone. Remove it.
        --nCoords;

        for (label ii = 0; ii < nCoords; ++ii)
        {
            coord[ii] = coord[ii+1];
            nDivs[ii] = nDivs[ii+1];
            ratio[ii] = ratio[ii+1];
            first[ii] = first[ii+1];
            last[ii]  = last[ii+1];
        }
        coord[0] = hard_min_;
    }

    width[0] = coord[1] - coord[0];
    first[0] = end_val(width[0],ratio[0],nDivs[0]);
    last[0] = end_val(width[0],1.0/ratio[0],nDivs[0]);


    // Now do the upper edge zone
    {
        const label np2 = nCoords - 2;

        ratio[np2] = pars.outerRatio;
        nDivs[np2] = pars.nEdgeLayers;

        coord[np2+1] =
            coord[np2] + last[np2-1] * rn1r1(pars.outerRatio, nDivs[np2]);

        width[np2] = coord[np2+1] - coord[np2];

        first[np2] = last[np2-1];
        last[np2]  = first[np2] * pow(pars.outerRatio, nDivs[np2] - 1);
    }


    // If we have a zone along the length of the plant that is much
    // longer than the width or height, then makePDRMeshBlocks does
    // not produce a very good outer boundary shape.
    // So here we divide the zone into several zones of roughly equal width.
    // A bit commplicated because of the increasing or decreasing step
    // sizes across the zone.

    if (max_zone > 0)
    {
        for (label ii = 0; ii < (nCoords-1); ++ii)
        {
            if (width[ii] > max_zone && nDivs[ii] > 1)
            {
                // No. of extra zones
                const label nExtra = width[ii] / max_zone - 1;

                // Make space for extra zones
                for (label ij = nCoords-1; ij > ii; --ij)
                {
                    width[ij+nExtra] = width[ij];
                    ratio[ij+nExtra] = ratio[ij];
                    first[ij+nExtra] = first[ij];
                    last[ij+nExtra]  = last[ij];
                    nDivs[ij+nExtra] = nDivs[ij];
                    coord[ij+nExtra] = coord[ij];
                }
                nCoords += nExtra;

                const scalar subWidth = width[ii] / (nExtra + 1);
                scalar bdy = coord[ii];

                last[ii+nExtra] = last[ii];
                nDivs[ii+nExtra] = nDivs[ii]; // will be decremented
                width[ii+nExtra] = width[ii]; // will be decremented

                // Current position
                scalar here = coord[ii];
                scalar step = first[ii];

                // Look for cell boundaries close to where we want
                // to put the zone boundaries
                for (label ij = ii; ij < ii+nExtra; ++ij)
                {
                    label ist = 0;
                    bdy += subWidth;
                    do
                    {
                        here += step;
                        step *= ratio[ii];
                        ++ist;
                    } while ((here + 0.5 * step) < bdy);

                    // Found a cell boundary at 'here' that is close to bdy

                    // Create this sub-zone
                    last[ij] = step / ratio[ii];
                    first[ij+1] = step;
                    coord[ij+1] = here;
                    width[ij] = coord[ij+1] - coord[ij];
                    nDivs[ij] = ist;
                    ratio[ij+1] = ratio[ij];

                    // Decrement what remains for the last sub-zone
                    nDivs[ii+nExtra] -= ist;
                    width[ii+nExtra] -= width[ij];
                }

                ii += nExtra;
            }
        }
    }


    if (verbose())
    {
        printf("Zone     Cell widths                  Ratios\n");
        printf("start    first  last           left   mid   right\n");
        for (label ii = 0; ii < (nCoords-1); ++ii)
        {
            printf
            (
                "%6.2f %6.2f %6.2f       %c %6.2f %6.2f %6.2f\n",
                coord[ii], first[ii], last[ii],
                marker[ii],
                first[ii]/last[ii-1], ratio[ii], first[ii+1]/last[ii]
            );
        }
        printf("%6.2f\n", coord[nCoords-1]);
    }


    // Copy back into a gridControl form

    PDRblock::gridControl grid;

    grid.resize(nCoords);
    for (label i = 0; i < nCoords; ++i)
    {
        grid[i] = coord[i];
    }
    for (label i = 0; i < nCoords-1; ++i)
    {
        grid.divisions_[i] = nDivs[i];
    }
    for (label i = 0; i < nCoords-1; ++i)
    {
        grid.expansion_[i] = ratio[i];
    }

    return grid;
}


// ************************************************************************* //
