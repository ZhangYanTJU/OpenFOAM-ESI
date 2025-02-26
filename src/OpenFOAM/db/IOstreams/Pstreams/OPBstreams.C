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

#include "UOPstream.H"
#include "OPstream.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::UOPBstream::UOPBstream
(
    DynamicList<char>& sendBuf,
    const int communicator,
    const bool sendAtDestruct,
    IOstreamOption::streamFormat fmt
)
:
    UOPstreamBase
    (
        UPstream::commsTypes::scheduled,    // irrelevant
        UPstream::masterNo(),   // irrelevant
        sendBuf,
        UPstream::msgType(),    // irrelevant
        communicator,
        sendAtDestruct,
        fmt
    )
{}


Foam::OPBstream::OPBstream
(
    const int communicator,
    IOstreamOption::streamFormat fmt
)
:
    Pstream(UPstream::commsTypes::scheduled),  // type is irrelevant
    UOPBstream
    (
        Pstream::transferBuf_,
        communicator,
        true,  // sendAtDestruct
        fmt
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::UOPBstream::~UOPBstream()
{
    if (sendAtDestruct_)
    {
        if (!bufferIPCsend())
        {
            FatalErrorInFunction
                << "Failed broadcast message of size " << sendBuf_.size()
                << Foam::abort(FatalError);
        }
    }
}


// ************************************************************************* //
