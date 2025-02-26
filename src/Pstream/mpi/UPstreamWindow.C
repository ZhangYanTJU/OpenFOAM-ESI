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

#include "PstreamGlobals.H"
#include "profilingPstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::UPstream::Window::Window() noexcept
:
    UPstream::Window(MPI_WIN_NULL)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::UPstream::Window::good() const noexcept
{
    return MPI_WIN_NULL != PstreamUtils::Cast::to_mpi(*this);
}


void Foam::UPstream::Window::reset() noexcept
{
    *this = UPstream::Window(MPI_WIN_NULL);
}


int Foam::UPstream::Window::size() const
{
    int val = 0;

    MPI_Win win = PstreamUtils::Cast::to_mpi(*this);
    MPI_Group group;

    // Get num of ranks from the group information
    if
    (
        (MPI_WIN_NULL != win)
     && (MPI_SUCCESS == MPI_Win_get_group(win, &group))
    )
    {
        if (MPI_SUCCESS != MPI_Group_size(group, &val))
        {
            val = 0;
        }
        MPI_Group_free(&group);
    }

    return val;
}


// ************************************************************************* //
