/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

Foam::UPstream::Request::Request() noexcept
:
    UPstream::Request(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::UPstream::Request::good() const noexcept
{
    return false;
}


void Foam::UPstream::Request::reset() noexcept
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::UPstream::nRequests() noexcept { return 0; }

void Foam::UPstream::resetRequests(const label n) {}

void Foam::UPstream::waitRequests(const label pos) {}
void Foam::UPstream::waitRequests(UList<UPstream::Request>&) {}

bool Foam::UPstream::waitAnyRequest(const label pos)
{
    return false;
}

bool Foam::UPstream::waitSomeRequests
(
    const label pos,
    DynamicList<int>* indices
)
{
    if (indices) indices->clear();
    return false;
}

Foam::label Foam::UPstream::waitAnyRequest(UList<UPstream::Request>&)
{
    return -1;
}

void Foam::UPstream::waitRequest(const label i) {}
void Foam::UPstream::waitRequest(UPstream::Request&) {}

bool Foam::UPstream::finishedRequest(const label i) { return true; }
bool Foam::UPstream::finishedRequest(UPstream::Request&) { return true; }

bool Foam::UPstream::finishedRequests(UList<UPstream::Request>&)
{
    return true;
}


bool Foam::UPstream::finishedRequestPair(label& req1, label& req2)
{
    req1 = -1;
    req2 = -1;
    return true;
}


void Foam::UPstream::waitRequestPair(label& req1, label& req2)
{
    req1 = -1;
    req2 = -1;
}


// ************************************************************************* //