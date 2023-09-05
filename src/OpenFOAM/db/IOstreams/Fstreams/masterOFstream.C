/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "masterOFstream.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "Pstream.H"
#include "masterUncollatedFileOperation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::masterOFstream::checkWrite
(
    const fileName& fName,
    const char* str,
    const std::streamsize len
)
{
    if (!str || !len)
    {
        // Can probably skip all of this if there is nothing to write
        return;
    }

    Foam::mkDir(fName.path());

    OFstream os
    (
        atomic_,
        fName,
        IOstreamOption(IOstreamOption::BINARY, version(), compression_),
        append_
    );
    if (!os.good())
    {
        FatalIOErrorInFunction(os)
            << "Could not open file " << fName << nl
            << exit(FatalIOError);
    }

    // Write characters directly to std::ostream
    os.writeRaw(str, len);

    if (!os.good())
    {
        FatalIOErrorInFunction(os)
            << "Failed writing to " << fName << nl
            << exit(FatalIOError);
    }
}


void Foam::masterOFstream::commit()
{
    // Take ownership of serialized content, without copying or reallocation
    DynamicList<char> charData(OCharStream::release());

    if (UPstream::parRun())
    {
        // Ignore content if not writing (reduces communication)
        if (!writeOnProc_)
        {
            charData.clear();
        }

        List<fileName> filePaths(UPstream::nProcs(comm_));
        filePaths[UPstream::myProcNo(comm_)] = pathName_;
        Pstream::gatherList(filePaths, UPstream::msgType(), comm_);

        // Test for identical output paths
        bool uniform =
        (
            UPstream::master(comm_)
          ? fileOperation::uniformFile(filePaths)
          : true
        );

        Pstream::broadcast(uniform, comm_);

        if (uniform)
        {
            // Identical file paths - write on master
            if (UPstream::master(comm_) && writeOnProc_)
            {
                checkWrite(pathName_, charData);
            }
            return;
        }

        // Different files
        // ---------------
        // Current strategy is to setup all non-blocking send/recv
        // using the probed message size to establish the recv size
        // (to avoid an additional communication of the sizes).
        //
        // For ranks with writeOnProc=false, the message size is 0.

        // An alternative approach would be to gather recv sizes
        // to avoid zero-sized messages and/or use double buffering
        // to recv into a buffer and write.
        //
        // const labelList recvSizes
        // (
        //     UPstream::listGatherValues<label>
        //     (
        //         (UPstream::is_subrank(comm_) ? charData.size() : label(0)),
        //         comm_
        //     )
        // );

        const label startOfRequests = UPstream::nRequests();

        // Some unique tag for this read/write/probe grouping
        const int messageTag = UPstream::msgType() + 256;

        if (UPstream::is_subrank(comm_))
        {
            // Send to master. When (!writeOnProc_) it is zero-sized.
            UOPstream::write
            (
                UPstream::commsTypes::nonBlocking,
                UPstream::masterNo(),
                charData.cdata_bytes(),
                charData.size_bytes(),
                messageTag,
                comm_
            );
        }
        else if (UPstream::master(comm_))
        {
            // The receive slots
            List<List<char>> procBuffers(UPstream::nProcs(comm_));

            const auto recvProcs = UPstream::subProcs(comm_);

            for (const int proci : recvProcs)
            {
                auto& procSlice = procBuffers[proci];

                // Probe the message size
                std::pair<int, int64_t> probed =
                    UPstream::probeMessage
                    (
                        UPstream::commsTypes::scheduled,  // blocking call
                        proci,
                        messageTag,
                        comm_
                    );

                procSlice.resize_nocopy(probed.second);

                // Receive content (can also be zero-sized)
                UIPstream::read
                (
                    UPstream::commsTypes::nonBlocking,
                    proci,
                    procSlice.data_bytes(),
                    procSlice.size_bytes(),
                    messageTag,
                    comm_
                );
            }

            if (writeOnProc_)
            {
                // Write non-empty master data
                checkWrite(pathName_, charData);
            }

            // Poll for completed receive requests and dispatch
            DynamicList<int> indices(recvProcs.size());
            while
            (
                UPstream::waitSomeRequests
                (
                    startOfRequests,
                    recvProcs.size(),
                   &indices
                )
            )
            {
                for (const int idx : indices)
                {
                    const int proci = recvProcs[idx];
                    auto& procSlice = procBuffers[proci];

                    if (!procSlice.empty())
                    {
                        // Write non-empty sub-proc data
                        checkWrite(filePaths[proci], procSlice);
                    }

                    // Eager cleanup?
                    // TBD: procSlice.clear();
                }
            }
        }

        UPstream::waitRequests(startOfRequests);
    }
    else
    {
        // Write (non-empty) data
        checkWrite(pathName_, charData);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::masterOFstream::masterOFstream
(
    IOstreamOption::atomicType atomic,
    const label comm,
    const fileName& pathName,
    IOstreamOption streamOpt,
    IOstreamOption::appendType append,
    const bool writeOnProc
)
:
    OCharStream(streamOpt),
    pathName_(pathName),
    atomic_(atomic),
    compression_(streamOpt.compression()),
    append_(append),
    writeOnProc_(writeOnProc),
    comm_(comm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::masterOFstream::~masterOFstream()
{
    commit();
}


// ************************************************************************* //
