/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2025 OpenCFD Ltd.
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

Description
    Gather data from all processors onto single processor according to some
    communication schedule (usually tree-to-master).
    The gathered data will be a single value constructed from the values
    on individual processors using a user-specified operator.

\*---------------------------------------------------------------------------*/

#include "IPstream.H"
#include "OPstream.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, class BinaryOp>
void Foam::Pstream::gather
(
    T& value,
    const BinaryOp& bop,
    const int tag,
    const label communicator
)
{
    if (UPstream::is_parallel(communicator))
    {
        // Communication order
        const auto& comms = UPstream::whichCommunication(communicator);
        // if (comms.empty()) return;  // extra safety?
        const auto& myComm = comms[UPstream::myProcNo(communicator)];

        // Receive from my downstairs neighbours
        for (const auto belowID : myComm.below())
        {
            T received;

            if constexpr (is_contiguous_v<T>)
            {
                UIPstream::read
                (
                    UPstream::commsTypes::scheduled,
                    belowID,
                    reinterpret_cast<char*>(&received),
                    sizeof(T),
                    tag,
                    communicator
                );
            }
            else
            {
                IPstream::recv(received, belowID, tag, communicator);
            }

            value = bop(value, received);
        }

        // Send up value
        if (myComm.above() >= 0)
        {
            if constexpr (is_contiguous_v<T>)
            {
                UOPstream::write
                (
                    UPstream::commsTypes::scheduled,
                    myComm.above(),
                    reinterpret_cast<const char*>(&value),
                    sizeof(T),
                    tag,
                    communicator
                );
            }
            else
            {
                OPstream::send(value, myComm.above(), tag, communicator);
            }
        }
    }
}


template<class T>
Foam::List<T> Foam::Pstream::listGatherValues
(
    const T& localValue,
    const label comm,
    [[maybe_unused]] const int tag
)
{
    if constexpr (is_contiguous_v<T>)
    {
        // UPstream version is contiguous only
        return UPstream::listGatherValues(localValue, comm);
    }
    else
    {
        List<T> allValues;

        if (UPstream::is_parallel(comm))
        {
            const label numProc = UPstream::nProcs(comm);

            if (UPstream::master(comm))
            {
                allValues.resize(numProc);

                // Non-trivial to manage non-blocking gather without a
                // PEX/NBX approach (eg, PstreamBuffers).
                // Leave with simple exchange for now

                allValues[0] = localValue;

                for (int proci = 1; proci < numProc; ++proci)
                {
                    IPstream::recv(allValues[proci], proci, tag, comm);
                }
            }
            else if (UPstream::is_rank(comm))
            {
                OPstream::send(localValue, UPstream::masterNo(), tag, comm);
            }
        }
        else
        {
            // non-parallel: return own value
            // TBD: only when UPstream::is_rank(comm) as well?
            allValues.resize(1);
            allValues[0] = localValue;
        }

        return allValues;
    }
}


template<class T>
T Foam::Pstream::listScatterValues
(
    const UList<T>& allValues,
    const label comm,
    [[maybe_unused]] const int tag
)
{
    if constexpr (is_contiguous_v<T>)
    {
        // UPstream version is contiguous only
        return UPstream::listScatterValues(allValues, comm);
    }
    else
    {
         T localValue{};

         if (UPstream::is_parallel(comm))
         {
             const label numProc = UPstream::nProcs(comm);

             if (UPstream::master(comm) && allValues.size() < numProc)
             {
                 FatalErrorInFunction
                     << "Attempting to send " << allValues.size()
                     << " values to " << numProc << " processors" << endl
                     << Foam::abort(FatalError);
             }

             if (UPstream::master(comm))
             {
                 const label startOfRequests = UPstream::nRequests();

                 List<DynamicList<char>> sendBuffers(numProc);

                 for (int proci = 1; proci < numProc; ++proci)
                 {
                     UOPstream toProc
                     (
                         UPstream::commsTypes::nonBlocking,
                         proci,
                         sendBuffers[proci],
                         tag,
                         comm
                     );
                     toProc << allValues[proci];
                 }

                 // Wait for outstanding requests
                 UPstream::waitRequests(startOfRequests);

                 return allValues[0];
             }
             else if (UPstream::is_rank(comm))
             {
                 IPstream::recv(localValue, UPstream::masterNo(), tag, comm);
             }
         }
         else
         {
             // non-parallel: return first value
             // TBD: only when UPstream::is_rank(comm) as well?

             if (!allValues.empty())
             {
                 return allValues[0];
             }
          }

          return localValue;
     }
}


// ************************************************************************* //
