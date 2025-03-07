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

#include "UPstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::UPstream::Window::Window() noexcept
:
    UPstream::Window(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::UPstream::Window::good() const noexcept
{
    return false;
}


void Foam::UPstream::Window::reset() noexcept
{}


int Foam::UPstream::Window::size() const
{
    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

std::pair<void*,int64_t>
Foam::UPstream::Window::mpi_win_allocate
(
    std::streamsize num_elements,
    int disp_unit,
    UPstream::Communicator communicator,
    const bool shared
)
{
    NotImplemented;
    return {nullptr, 0};
}


std::pair<void*,int64_t>
Foam::UPstream::Window::mpi_win_allocate
(
    std::streamsize num_elements,
    int disp_unit,
    int communicator,
    const bool shared
)
{
    NotImplemented;
    return {nullptr, 0};
}


bool Foam::UPstream::Window::mpi_win_create
(
    void *baseptr,
    std::streamsize num_elements,
    int disp_unit,
    UPstream::Communicator communicator
)
{
    NotImplemented;
    return false;
}


bool Foam::UPstream::Window::mpi_win_create
(
    void *baseptr,
    std::streamsize num_elements,
    int disp_unit,
    int communicator
)
{
    NotImplemented;
    return false;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

void Foam::UPstream::Window::close()
{}


// * * * * * * * * * * * * * * * Synchronization * * * * * * * * * * * * * * //

void Foam::UPstream::Window::mpi_win_flushing(int rank, bool local)
{}


void Foam::UPstream::Window::sync()
{}


void Foam::UPstream::Window::mpi_win_locking(int rank, bool exclusive)
{}


void Foam::UPstream::Window::mpi_win_unlocking(int rank)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::UPstream::Window::get_data
(
    void* origin,
    std::streamsize count,
    const UPstream::dataTypes dataTypeId,
    int target_rank,
    int target_disp
) const
{
    NotImplemented;
    return false;
}


bool Foam::UPstream::Window::put_data
(
    const void* origin,
    std::streamsize count,
    const UPstream::dataTypes dataTypeId,
    int target_rank,
    int target_disp
) const
{
    NotImplemented;
    return false;
}


bool Foam::UPstream::Window::put_data
(
    const UPstream::opCodes opCodeId,
    const void* origin,
    std::streamsize count,
    const UPstream::dataTypes dataTypeId,
    int target_rank,
    int target_disp
) const
{
    NotImplemented;
    return false;
}


bool Foam::UPstream::Window::mpi_fetch_and_op
(
    const UPstream::opCodes opCodeId,
    const void* origin,
    void* result,
    const UPstream::dataTypes dataTypeId,
    int target_rank,
    int target_disp
) const
{
    NotImplemented;
    return false;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::UPstream::Window::is_shared(const bool failNonShared) const
{
    return false;
}


std::pair<void*,int64_t>
Foam::UPstream::Window::mpi_win_query
(
    UPstream::Window window,
    const int expected_disp_unit
)
{
    // No window to query
    NotImplemented;
    return {nullptr, 0};
}


std::pair<void*,int64_t>
Foam::UPstream::Window::mpi_win_query_shared
(
    UPstream::Window window,
    int target_rank,
    const int expected_disp_unit
)
{
    // No window to query
    NotImplemented;
    return {nullptr, 0};
}


// ************************************************************************* //
