/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2021-2025 OpenCFD Ltd.
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
    Foam::OPstream

Description
    Output inter-processor communications stream.

SourceFiles
    OPstreams.C

\*---------------------------------------------------------------------------*/

#include "Pstream.H"

#ifndef Foam_OPstream_H
#define Foam_OPstream_H

#include "UOPstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class OPstream Declaration
\*---------------------------------------------------------------------------*/

class OPstream
:
    public Pstream,
    public UOPstream
{
public:

    // Constructors

        //- Construct for given process index to send to
        OPstream
        (
            const UPstream::commsTypes commsType,
            const int toProcNo,
            const int bufferSize = 0,  //!< optional buffer size
            const int tag = UPstream::msgType(),
            const int communicator = UPstream::worldComm,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        );


    // Static Functions

        //- Serialize a value and send (buffered/blocking or standard mode).
        //- Uses \c operator<< for serialization
        template<class Type>
        static void send
        (
            const Type& value,
            //! buffered or scheduled only!
            const UPstream::commsTypes commsType,
            const int toProcNo,
            const int tag = UPstream::msgType(),
            const int communicator = UPstream::worldComm,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        )
        {
            OPstream os(commsType, toProcNo, 0, tag, communicator, fmt);
            os << value;
        }

        //- Serialize a value and send (standard mode).
        //- Uses \c operator<< for serialization
        template<class Type>
        static void send
        (
            const Type& value,
            const int toProcNo,
            const int tag = UPstream::msgType(),
            const int communicator = UPstream::worldComm,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        )
        {
            OPstream::send
            (
                value,
                UPstream::commsTypes::scheduled,  // ie, MPI_Send()
                toProcNo,
                tag,
                communicator,
                fmt
            );
        }
};


/*---------------------------------------------------------------------------*\
                        Class OPBstream Declaration
\*---------------------------------------------------------------------------*/

//- Output inter-processor communications stream
//- using MPI broadcast.
class OPBstream
:
    public Pstream,
    public UOPBstream
{
public:

    // Constructors

        //- Construct with optional communicator and write format.
        explicit OPBstream
        (
            const int communicator = UPstream::worldComm,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        );


    // Static Functions

        //- Use all send methods from base
        using UOPBstream::send;

        //- Serialize a value and broadcast (root == UPstream::masterNo()).
        //- Uses \c operator<< for serialization
        template<class Type>
        static void send
        (
            const Type& value,
            const int communicator = UPstream::worldComm
        )
        {
            OPBstream os(communicator);
            {
                os << value;
            }
        }

        //- Serialize multiple items and broadcast the buffer
        //- Uses \c operator<< for serialization
        template<class Type, class... Args>
        static void sends
        (
            const int communicator,
            Type& value,
            Args&&... values
        )
        {
            OPBstream os(communicator);
            {
                Detail::outputLoop(os, value, std::forward<Args>(values)...);
                // Depending on compiler support:
                // Pack via fold expression
                // (((os << value) << std::forward<Args>(values)), ...);
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
