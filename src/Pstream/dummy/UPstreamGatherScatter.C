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

#include "UPstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::UPstream::mpi_gather
(
    const void* sendData,
    void* recvData,
    int count,
    const UPstream::dataTypes dataTypeId,

    const int communicator,
    UPstream::Request* req
)
{}


void Foam::UPstream::mpi_scatter
(
    const void* sendData,
    void* recvData,
    int count,
    const UPstream::dataTypes dataTypeId,

    const int communicator,
    UPstream::Request* req
)
{}


void Foam::UPstream::mpi_allgather
(
    void* allData,
    int count,
    const UPstream::dataTypes dataTypeId,

    const int communicator,
    UPstream::Request* req
)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::UPstream::mpi_gatherv
(
    const void* sendData,
    int sendCount,
    void* recvData,
    const UList<int>& recvCounts,
    const UList<int>& recvOffsets,

    const UPstream::dataTypes dataTypeId,
    const int communicator
)
{}


void Foam::UPstream::mpi_scatterv
(
    const void* sendData,
    const UList<int>& sendCounts,
    const UList<int>& sendOffsets,

    void* recvData,
    int recvCount,

    const UPstream::dataTypes dataTypeId,
    const int communicator
)
{}


// ************************************************************************* //
