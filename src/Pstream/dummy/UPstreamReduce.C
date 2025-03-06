/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2025 OpenCFD Ltd.
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

#include "Pstream.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// Special reductions for bool

void Foam::UPstream::reduceAnd(bool& value, const int communicator)
{}

void Foam::UPstream::reduceOr(bool& value, const int communicator)
{}


void Foam::reduce
(
    bool& value,
    Foam::andOp<bool>,
    const int tag,
    const int communicator
)
{}

void Foam::reduce
(
    bool& value,
    Foam::orOp<bool>,
    const int tag,
    const int communicator
)
{}


// * * * * * * * * * * Protected Static Member Functions * * * * * * * * * * //

void Foam::UPstream::mpi_reduce
(
    void* values,
    int count,
    const UPstream::dataTypes dataTypeId,
    const UPstream::opCodes opCodeId,
    const int communicator,
    UPstream::Request* req
)
{}


void Foam::UPstream::mpi_allreduce
(
    void* values,
    int count,
    const UPstream::dataTypes dataTypeId,
    const UPstream::opCodes opCodeId,
    const int communicator,
    UPstream::Request* req
)
{}


// ************************************************************************* //
