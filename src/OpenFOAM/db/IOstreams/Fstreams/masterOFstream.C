/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
#include "PstreamBuffers.H"
#include "masterUncollatedFileOperation.H"
#include <algorithm>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::masterOFstream::checkWrite
(
    const fileName& fName,
    const char* str,
    std::streamsize len
)
{
    if (!len)
    {
        // Can probably skip all of this if there is nothing to write
        return;
    }

    mkDir(fName.path());

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

    // Use writeRaw() instead of writeQuoted(string,false) to output
    // characters directly.

    os.writeRaw(str, len);

    if (!os.good())
    {
        FatalIOErrorInFunction(os)
            << "Failed writing to " << fName << nl
            << exit(FatalIOError);
    }
}


void Foam::masterOFstream::checkWrite
(
    const fileName& fName,
    const std::string& s
)
{
    checkWrite(fName, s.data(), s.length());
}


void Foam::masterOFstream::commit()
{
    if (Pstream::parRun())
    {
        List<fileName> filePaths(Pstream::nProcs(comm_));
        filePaths[Pstream::myProcNo(comm_)] = pathName_;
        Pstream::gatherList(filePaths, UPstream::msgType(), comm_);

        bool uniform =
        (
            Pstream::master(comm_) && fileOperation::uniformFile(filePaths)
        );

        Pstream::broadcast(uniform, comm_);

        if (uniform)
        {
            if (Pstream::master(comm_) && valid_)
            {
                checkWrite(pathName_, this->str());
            }

            this->reset();
            return;
        }

        boolList procValid(UPstream::listGatherValues<bool>(valid_, comm_));

        // Different files
        PstreamBuffers pBufs(comm_, Pstream::commsTypes::nonBlocking);

        // Send my buffer to master
        if (!Pstream::master(comm_))
        {
            UOPstream os(Pstream::masterNo(), pBufs);
            string s(this->str());
            this->reset();

            os.write(s.data(), s.length());
        }

        labelList recvSizes;
        pBufs.finishedGathers(recvSizes);

        if (Pstream::master(comm_))
        {
            // Write master data
            if (procValid[Pstream::masterNo()])
            {
                checkWrite(filePaths[Pstream::masterNo()], this->str());
            }
            this->reset();

            // Find the max receive size
            recvSizes[Pstream::masterNo()] = 0;
            List<char> buf
            (
                *std::max_element(recvSizes.cbegin(), recvSizes.cend())
            );

            for (const int proci : Pstream::subProcs(comm_))
            {
                UIPstream is(proci, pBufs);

                const std::streamsize count(recvSizes[proci]);
                is.read(buf.data(), count);

                if (procValid[proci])
                {
                    checkWrite(filePaths[proci], buf.cdata(), count);
                }
            }
        }
    }
    else
    {
        checkWrite(pathName_, this->str());
        this->reset();
    }

    // This method is only called once (internally)
    // so no need to clear/flush old buffered data
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::masterOFstream::masterOFstream
(
    IOstreamOption::atomicType atomic,
    const label comm,
    const fileName& pathName,
    IOstreamOption streamOpt,
    IOstreamOption::appendType append,
    const bool valid
)
:
    OStringStream(streamOpt),
    pathName_(pathName),
    atomic_(atomic),
    compression_(streamOpt.compression()),
    append_(append),
    valid_(valid),
    comm_(comm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::masterOFstream::~masterOFstream()
{
    commit();
}


// ************************************************************************* //
