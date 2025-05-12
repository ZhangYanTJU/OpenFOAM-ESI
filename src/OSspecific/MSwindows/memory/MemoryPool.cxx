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

#include "MemoryPool.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

// bool Foam::MemoryPool::create(const std::string& ctrl, bool verbose)
// {
//     return false;
// }


bool Foam::MemoryPool::create(bool verbose)
{
    // No banner information since it is currently never an option
    return false;
}


void Foam::MemoryPool::destroy(bool verbose)
{}


bool Foam::MemoryPool::active() noexcept
{
    return false;
}


bool Foam::MemoryPool::suspend() noexcept
{
    return false;
}


void Foam::MemoryPool::resume() noexcept
{}


bool Foam::MemoryPool::is_pool(void* ptr)
{
    return false;
}


void* Foam::MemoryPool::try_allocate(std::size_t nbytes)
{
    return nullptr;
}


bool Foam::MemoryPool::try_deallocate(void* ptr)
{
    return (!ptr);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
