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

#include "PDRfitMeshScans.H"
#include "PDRobstacle.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::PDRfitMeshScans::prepare
(
    const UList<PDRobstacle>& obstacles,
    const PDRfitMeshParams& fitParams,
    scalar cellWidth
)
{
    reset();

    // Scan for obstacle limits
    scanLimits(obstacles);

    // Clip with optional user bounds

    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        scalarMinMax& limits = (*this)[cmpt].limits();

        if (fitParams.minBounds[cmpt].has_value())
        {
            limits.min() = fitParams.minBounds[cmpt].value();
        }
        if (fitParams.maxBounds[cmpt].has_value())
        {
            limits.max() = fitParams.maxBounds[cmpt].value();
        }
    }

    // Get obstacle locations and subgrid length (if any)
    scanAreas(obstacles, fitParams, cellWidth);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::PDRfitMeshScans::volume() const
{
    scalar vol = 1;

    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        const scalarMinMax& limits = (*this)[cmpt].limits();

        if (limits.valid())
        {
            vol *= limits.mag();
        }
        else
        {
            return 0;
        }
    }

    return vol;
}


Foam::scalar Foam::PDRfitMeshScans::cellVolume() const
{
    scalar vol = 1;

    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        vol *= (*this)[cmpt].stepWidth();
    }

    return vol;
}


void Foam::PDRfitMeshScans::reset()
{
    subgridLen_ = 0;

    for (PDRfitMeshScan& pln : *this)
    {
        pln.reset();
    }
}


void Foam::PDRfitMeshScans::resize
(
    const scalar cellWidth,
    const PDRfitMeshParams& pars
)
{
    for (PDRfitMeshScan& pln : *this)
    {
        pln.resize(cellWidth, pars);
    }
}




void Foam::PDRfitMeshScans::print(Ostream& os) const
{
    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        os << vector::componentNames[cmpt] << ' ';
        (*this)[cmpt].print(os);
    }
}


void Foam::PDRfitMeshScans::scanLimits
(
    const UList<PDRobstacle>& obstacles
)
{
    for (const PDRobstacle& obs : obstacles)
    {
        switch (obs.typeId)
        {
            case PDRobstacle::CYLINDER:
            case PDRobstacle::DIAG_BEAM:
            {
                (*this)[obs.orient].adjustLimits
                (
                    obs.pt[obs.orient],
                    obs.pt[obs.orient] + obs.len()
                );
                break;
            }

            case PDRobstacle::CUBOID_1:
            case PDRobstacle::LOUVRE_BLOWOFF:
            case PDRobstacle::CUBOID:
            case PDRobstacle::WALL_BEAM:
            case PDRobstacle::GRATING:
            case PDRobstacle::RECT_PATCH:
            {
                (*this).x().adjustLimits(obs.x(), obs.x() + obs.span.x());
                (*this).y().adjustLimits(obs.y(), obs.y() + obs.span.y());
                (*this).z().adjustLimits(obs.z(), obs.z() + obs.span.z());
                break;
            }
        }
    }
}


Foam::scalar
Foam::PDRfitMeshScans::scanAreas
(
    const UList<PDRobstacle>& obstacles,
    const scalar minSubgridLen
)
{
    scalar subgridLen = 0;

    const bool checkSubgrid = (minSubgridLen > 0);

    if (verbose())
    {
        Info<< "Scan " << obstacles.size()
            << " obstacles. Check for subgrid < " << minSubgridLen
            << nl;
    }

    // When sorting flat elements, convert vector -> List
    List<scalar> gridSpan(3);

    for (const PDRobstacle& obs : obstacles)
    {
        switch (obs.typeId)
        {
            case PDRobstacle::CYLINDER:
            {
                // #ifdef FULLDEBUG
                // Info<< "Add " << obs.info() << nl;
                // #endif

                (*this)[obs.orient].addArea
                (
                    sqr(0.5*obs.dia()),
                    obs.pt[obs.orient],
                    obs.pt[obs.orient] + obs.len()
                );

                // Store total length of subgrid obstcales
                if (checkSubgrid && obs.dia() < minSubgridLen)
                {
                    subgridLen += obs.len();
                }
                break;
            }

            case PDRobstacle::DIAG_BEAM:
            {
                // #ifdef FULLDEBUG
                // Info<< "Add " << obs.info() << nl;
                // #endif

                (*this)[obs.orient].addArea
                (
                    (obs.wa * obs.wb),
                    obs.pt[obs.orient],
                    obs.pt[obs.orient] + obs.len()
                );

                // Store total length of subgrid obstcales
                if (checkSubgrid && Foam::max(obs.wa, obs.wb) < minSubgridLen)
                {
                    subgridLen += obs.len();
                }
                break;
            }

            case PDRobstacle::CUBOID_1:
            case PDRobstacle::LOUVRE_BLOWOFF:
            case PDRobstacle::CUBOID:
            case PDRobstacle::WALL_BEAM:
            case PDRobstacle::GRATING:
            case PDRobstacle::RECT_PATCH:
            {
                // #ifdef FULLDEBUG
                // Info<< "Add " << obs.info() << nl;
                // #endif

                // Can be lazy here, addArea() filters out zero areas

                (*this)[vector::X].addArea
                (
                    (obs.span.y() * obs.span.z()),
                    obs.x(),
                    obs.x() + obs.span.x()
                );

                (*this)[vector::Y].addArea
                (
                    (obs.span.z() * obs.span.x()),
                    obs.y(),
                    obs.y() + obs.span.y()
                );

                (*this)[vector::Z].addArea
                (
                    (obs.span.x() * obs.span.y()),
                    obs.z(),
                    obs.z() + obs.span.z()
                );

                if (checkSubgrid)
                {
                    gridSpan[0] = obs.span.x();
                    gridSpan[1] = obs.span.y();
                    gridSpan[2] = obs.span.z();
                    Foam::sort(gridSpan);

                    // Ignore zero dimension when considering subgrid

                    // Store total length of subgrid obstcales
                    if (gridSpan[1] < minSubgridLen)
                    {
                        subgridLen += gridSpan.last();
                    }
                }
                break;
            }
        }
    }

    return subgridLen;
}


void Foam::PDRfitMeshScans::scanAreas
(
    const UList<PDRobstacle>& obstacles
)
{
    (void)scanAreas(obstacles, -GREAT);
}


void Foam::PDRfitMeshScans::scanAreas
(
    const UList<PDRobstacle>& obstacles,
    const PDRfitMeshParams& fitParams,
    scalar cellWidth
)
{
    resize(mag(cellWidth), fitParams);

    const scalar minSubgridLen =
        (cellWidth < 0 ? (mag(cellWidth) * fitParams.widthFactor) : 0);

    subgridLen_ = scanAreas(obstacles, minSubgridLen);
}


Foam::Vector<Foam::PDRblock::gridControl>
Foam::PDRfitMeshScans::calcGriding
(
    const UList<PDRobstacle>& obstacles,
    const PDRfitMeshParams& fitParams,
    scalar cellWidth
)
{
    prepare(obstacles, fitParams, cellWidth);

    const scalar innerVol = volume();

    // Set hard lower limit at the ground

    z().hard_min(fitParams.ground);

    if (z().limits().min() < fitParams.ground)
    {
        z().limits().min() = fitParams.ground;
    }

    if (verbose())
    {
        print(Info);
        Info<< "cellWidth = " << cellWidth << nl;
        Info<< "inner volume: " << innerVol << nl;
        Info<< "subgrid length: " << subgridLen_ << nl;
    }


    // Read in obstacles and record face planes Also record length of
    // subgrid obstacles so that we can optimise cell size to get
    // desired average no, of obstacles oper cell. Determinuing which
    // obstacles are subgrid is iself dependent on the cell size.
    // Therefore we iterate (if cell_width is -ve) If cell_width is
    // +ve just use the user-supplied value.

    {
        // The cell-widths during optimisation
        scalar prev_cw = 0, cw = 0;

        const bool optimiseWidth = (cellWidth < 0);

        for (label nIter = 0; nIter < fitParams.maxIterations; ++nIter)
        {
            scanAreas(obstacles, fitParams, cellWidth);

            if (cellWidth < 0)
            {
                constexpr scalar relax = 0.7;

                prev_cw = cw;

                if (subgridLen_ < cellWidth)
                {
                    FatalErrorInFunction
                        << "No sub-grid obstacles found" << endl
                        << exit(FatalError);
                }

                // Optimise to average subgrid obstacles per cell
                cw = sqrt(fitParams.obsPerCell * innerVol / subgridLen_);
                cw = Foam::min(cw, fitParams.maxCellWidth);

                if (cw > cellWidth * fitParams.maxWidthEstimate)
                {
                    Warning
                        << "Calculated cell width more than"
                           " maxWidthEstimate x estimate." << nl
                        << "Too few sub-grid obstacles?"
                        << nl;
                }

                Info<< "Current cellwidth: " << cw << nl;

                scalar ratio = cw / cellWidth;

                if (ratio < 1)
                {
                    ratio = 1/ratio;
                }

                if (ratio < fitParams.maxWidthRatio)
                {
                    cellWidth = -cellWidth;
                }
                else
                {
                    cellWidth = relax * (-cw) + (1.0 - relax) * cellWidth;
                }
                Info<< "Cell width: " << cellWidth << nl;
            }

            if (cellWidth > 0)
            {
                break;
            }
        }

        if (cellWidth < 0)
        {
            cellWidth = 0.5 * (cw + prev_cw);
        }

        if (optimiseWidth)
        {
            Info<< nl << "Final cell width: " << cellWidth << nl << nl;
        }
    }


    // Now we fit the mesh to the planes and write out the result

    resize(cellWidth, fitParams);
    cellWidth = cbrt(cellVolume()) / fitParams.areaWidthFactor;

    scanAreas(obstacles, fitParams, cellWidth);

    const scalar max_zone =
    (
        (z().limits().mag() + fitParams.nEdgeLayers * cellWidth)
      * fitParams.maxZoneToHeight
    );


    Vector<PDRblock::gridControl> griding;

    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        griding[cmpt] =
            (*this)[cmpt].calcGridControl(cellWidth, max_zone, fitParams);
    }

    return griding;
}


// ************************************************************************* //
