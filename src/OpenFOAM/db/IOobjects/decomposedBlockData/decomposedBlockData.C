/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "decomposedBlockData.H"
#include "OPstream.H"
#include "IPstream.H"
#include "PstreamBuffers.H"
#include "Fstream.H"
#include "dictionary.H"
#include "objectRegistry.H"
#include "SubList.H"
#include "charList.H"
#include "labelPair.H"
#include "masterUncollatedFileOperation.H"
#include "ListStream.H"
#include "StringStream.H"
#include "registerSwitch.H"
#include <cctype>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(decomposedBlockData, 0);
}


int Foam::decomposedBlockData::collatedBlockFormat_
(
    Foam::debug::optimisationSwitch("collatedFormat", 1)
);

registerOptSwitch
(
    "collatedFormat",
    int,
    Foam::decomposedBlockData::collatedBlockFormat_
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decomposedBlockData::decomposedBlockData
(
    const label comm,
    const IOobject& io,
    const UPstream::commsTypes commsType
)
:
    regIOobject(io),
    commsType_(commsType),
    comm_(comm),
    streamOption_(),
    contentData_()
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningInFunction
            << "decomposedBlockData " << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but decomposedBlockData does not support automatic rereading."
            << endl;
    }
    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        read();
    }
}


// * * * * * * * * * * * * * * * Members Functions * * * * * * * * * * * * * //

bool Foam::decomposedBlockData::readBlockEntry
(
    Istream& is,
    List<char>& charData
)
{
    // Handle any of these:

    // 0.  NCHARS (...)
    // 1.  List<char> NCHARS (...)
    // 2.  processorN  List<char> NCHARS (...) ;

    is.fatalCheck(FUNCTION_NAME);
    token tok(is);
    is.fatalCheck(FUNCTION_NAME);

    // Dictionary format has primitiveEntry keyword:
    const bool isDictFormat = (tok.isWord() && !tok.isCompound());

    if (!isDictFormat && tok.good())
    {
        is.putBack(tok);
    }
    charData.readList(is);

    if (isDictFormat)
    {
        is.fatalCheck(FUNCTION_NAME);
        is >> tok;
        is.fatalCheck(FUNCTION_NAME);

        // Swallow trailing ';'
        if (tok.good() && !tok.isPunctuation(token::END_STATEMENT))
        {
            is.putBack(tok);
        }
    }

    return true;
}


std::streamoff Foam::decomposedBlockData::writeBlockEntry
(
    OSstream& os,
    const label blocki,
    const UList<char>& charData
)
{
    // Offset to the beginning of this output

    std::streamoff blockOffset = os.stdStream().tellp();

    const word procName("processor" + Foam::name(blocki));

    if (decomposedBlockData::oldBlockFormat())
    {
        // Bare old style
        os  << nl << "// " << procName << nl;
        charData.writeList(os) << nl;
    }
    else
    {
        // Dictionary format
        os  << nl;
        charData.writeEntry(procName, os);
    }

    return blockOffset;
}


std::streamoff Foam::decomposedBlockData::writeBlockEntry
(
    OSstream& os,
    IOstreamOption streamOptData,
    const regIOobject& io,
    const label blocki,
    const bool withLocalHeader
)
{
    // String(s) from all data to write
    string contentChars;
    {
        OStringStream os(streamOptData);

        bool ok = true;

        // Generate FoamFile header on master, without comment banner
        if (withLocalHeader)
        {
            const bool old = IOobject::bannerEnabled(false);

            ok = io.writeHeader(os);

            IOobject::bannerEnabled(old);
        }

        // Write the data to the Ostream
        ok = ok && io.writeData(os);

        if (!ok)
        {
            return std::streamoff(-1);
        }

        contentChars = os.str();
    }

    // The character data
    UList<char> charData
    (
        const_cast<char*>(contentChars.data()),
        label(contentChars.size())
    );

    return decomposedBlockData::writeBlockEntry(os, blocki, charData);
}


Foam::autoPtr<Foam::ISstream>
Foam::decomposedBlockData::readBlock
(
    const label blocki,
    ISstream& is,
    IOobject& headerIO
)
{
    if (debug)
    {
        Pout<< "decomposedBlockData::readBlock:"
            << " stream:" << is.name() << " attempt to read block " << blocki
            << endl;
    }

    // Extracted header information
    IOstreamOption streamOptData;
    unsigned streamLabelByteSize = sizeof(label);
    unsigned streamScalarByteSize = sizeof(scalar);

    is.fatalCheck(FUNCTION_NAME);
    token tok(is);
    is.fatalCheck(FUNCTION_NAME);

    const bool hasBlockHeader = tok.isWord(decomposedBlockData::headerName);

    if (hasBlockHeader)
    {
        // FoamBlock: header for block data
        dictionary headerDict(is, false);  // Read sub-dictionary content

        streamOptData = headerIO.parseHeader(headerDict);
        streamLabelByteSize = headerIO.labelByteSize();
        streamScalarByteSize = headerIO.scalarByteSize();
    }
    else
    {
        is.putBack(tok);
    }

    autoPtr<ISstream> realIsPtr;

    // Read master for header
    List<char> data;
    decomposedBlockData::readBlockEntry(is, data);

    if (blocki == 0)
    {
        realIsPtr.reset(new IListStream(std::move(data)));
        realIsPtr->name() = is.name();

        if (!hasBlockHeader)
        {
            // No separate FoamBlock:
            // Read header from first block,
            // advancing the stream position
            if (!headerIO.readHeader(*realIsPtr))
            {
                FatalIOErrorInFunction(*realIsPtr)
                    << "Problem while reading header for object "
                    << is.name() << nl
                    << exit(FatalIOError);
            }
        }
        else
        {
            // Apply stream settings
            realIsPtr().format(streamOptData.format());
            realIsPtr().version(streamOptData.version());
            realIsPtr().setLabelByteSize(streamLabelByteSize);
            realIsPtr().setScalarByteSize(streamScalarByteSize);
        }
    }
    else
    {
        if (!hasBlockHeader)
        {
            // No separate FoamBlock:
            // Read header from first block
            UIListStream headerStream(data);
            if (!headerIO.readHeader(headerStream))
            {
                FatalIOErrorInFunction(headerStream)
                    << "Problem while reading header for object "
                    << is.name() << nl
                    << exit(FatalIOError);
            }
            streamOptData = static_cast<IOstreamOption>(headerStream);
            streamLabelByteSize = headerStream.labelByteSize();
            streamScalarByteSize = headerStream.scalarByteSize();
        }

        for (label i = 1; i < blocki+1; i++)
        {
            // Read and discard data, only retain the last one
            decomposedBlockData::readBlockEntry(is, data);
        }
        realIsPtr.reset(new IListStream(std::move(data)));
        realIsPtr->name() = is.name();

        // Apply stream settings
        realIsPtr().format(streamOptData.format());
        realIsPtr().version(streamOptData.version());
        realIsPtr().setLabelByteSize(streamLabelByteSize);
        realIsPtr().setScalarByteSize(streamScalarByteSize);
    }

    return realIsPtr;
}


bool Foam::decomposedBlockData::readBlocks
(
    const label comm,
    autoPtr<ISstream>& isPtr,
    List<char>& data,
    const UPstream::commsTypes commsType
)
{
    if (debug)
    {
        Pout<< "decomposedBlockData::readBlocks:"
            << " stream:" << (isPtr ? isPtr->name() : "invalid")
            << " commsType:" << Pstream::commsTypeNames[commsType]
            << " comm:" << comm << endl;
    }

    bool ok = false;

    if (UPstream::master(comm))
    {
        auto& is = *isPtr;
        is.fatalCheck(FUNCTION_NAME);

        // Read master data
        decomposedBlockData::readBlockEntry(is, data);
    }

    if (commsType == UPstream::commsTypes::scheduled)
    {
        if (UPstream::master(comm))
        {
            // Master data already read ...
            auto& is = *isPtr;
            is.fatalCheck(FUNCTION_NAME);

            // Read and transmit slave data
            for (const int proci : UPstream::subProcs(comm))
            {
                List<char> elems;
                decomposedBlockData::readBlockEntry(is, elems);

                OPstream os
                (
                    UPstream::commsTypes::scheduled,
                    proci,
                    0,
                    UPstream::msgType(),
                    comm
                );
                os << elems;
            }

            ok = is.good();
        }
        else
        {
            IPstream is
            (
                UPstream::commsTypes::scheduled,
                UPstream::masterNo(),
                0,
                UPstream::msgType(),
                comm
            );
            is >> data;
        }
    }
    else
    {
        PstreamBuffers pBufs
        (
            UPstream::commsTypes::nonBlocking,
            UPstream::msgType(),
            comm
        );

        if (UPstream::master(comm))
        {
            // Master data already read ...
            auto& is = *isPtr;
            is.fatalCheck(FUNCTION_NAME);

            // Read and transmit slave data
            for (const int proci : UPstream::subProcs(comm))
            {
                List<char> elems;
                decomposedBlockData::readBlockEntry(is, elems);

                UOPstream os(proci, pBufs);
                os << elems;
            }
        }

        labelList recvSizes;
        pBufs.finishedSends(recvSizes);

        if (!UPstream::master(comm))
        {
            UIPstream is(UPstream::masterNo(), pBufs);
            is >> data;
        }
    }

    Pstream::scatter(ok, Pstream::msgType(), comm);

    return ok;
}


Foam::autoPtr<Foam::ISstream> Foam::decomposedBlockData::readBlocks
(
    const label comm,
    const fileName& fName,
    autoPtr<ISstream>& isPtr,
    IOobject& headerIO,
    const UPstream::commsTypes commsType
)
{
    if (debug)
    {
        Pout<< "decomposedBlockData::readBlocks:"
            << " stream:" << (isPtr ? isPtr->name() : "invalid")
            << " commsType:" << Pstream::commsTypeNames[commsType] << endl;
    }

    bool ok = false;
    List<char> data;
    autoPtr<ISstream> realIsPtr;

    if (UPstream::master(comm))
    {
        auto& is = *isPtr;

        // Extracted header information
        IOstreamOption streamOptData;
        unsigned streamLabelByteSize = sizeof(label);
        unsigned streamScalarByteSize = sizeof(scalar);

        is.fatalCheck(FUNCTION_NAME);
        token tok(is);
        is.fatalCheck(FUNCTION_NAME);

        const bool hasBlockHeader = tok.isWord(decomposedBlockData::headerName);

        if (hasBlockHeader)
        {
            // FoamBlock: header for block data
            dictionary headerDict(is, false);  // Read sub-dictionary content

            streamOptData = headerIO.parseHeader(headerDict);
            streamLabelByteSize = headerIO.labelByteSize();
            streamScalarByteSize = headerIO.scalarByteSize();
        }
        else
        {
            is.putBack(tok);
        }

        // Read master data
        decomposedBlockData::readBlockEntry(is, data);

        realIsPtr.reset(new IListStream(std::move(data)));
        realIsPtr->name() = fName;

        if (!hasBlockHeader)
        {
            // No separate FoamBlock:
            // Read header from first block,
            // advancing the stream position
            if (!headerIO.readHeader(*realIsPtr))
            {
                FatalIOErrorInFunction(*realIsPtr)
                    << "Problem while reading header for object "
                    << is.name() << nl
                    << exit(FatalIOError);
            }
        }
        else
        {
            // Apply stream settings
            realIsPtr().format(streamOptData.format());
            realIsPtr().version(streamOptData.version());
            realIsPtr().setLabelByteSize(streamLabelByteSize);
            realIsPtr().setScalarByteSize(streamScalarByteSize);
        }
    }

    if (commsType == UPstream::commsTypes::scheduled)
    {
        if (UPstream::master(comm))
        {
            // Master data already read ...
            auto& is = *isPtr;
            is.fatalCheck(FUNCTION_NAME);

            // Read and transmit slave data
            for (const int proci : UPstream::subProcs(comm))
            {
                decomposedBlockData::readBlockEntry(is, data);

                OPstream os
                (
                    UPstream::commsTypes::scheduled,
                    proci,
                    0,
                    UPstream::msgType(),
                    comm
                );
                os << data;
            }

            ok = is.good();
        }
        else
        {
            IPstream is
            (
                UPstream::commsTypes::scheduled,
                UPstream::masterNo(),
                0,
                UPstream::msgType(),
                comm
            );
            is >> data;

            realIsPtr.reset(new IListStream(std::move(data)));
            realIsPtr->name() = fName;
        }
    }
    else
    {
        PstreamBuffers pBufs
        (
            UPstream::commsTypes::nonBlocking,
            UPstream::msgType(),
            comm
        );

        if (UPstream::master(comm))
        {
            // Master data already read ...
            auto& is = *isPtr;
            is.fatalCheck(FUNCTION_NAME);

            // Read and transmit slave data
            for (const int proci : UPstream::subProcs(comm))
            {
                List<char> elems;
                decomposedBlockData::readBlockEntry(is, elems);

                UOPstream os(proci, pBufs);
                os << elems;
            }

            ok = is.good();
        }

        labelList recvSizes;
        pBufs.finishedSends(recvSizes);

        if (!UPstream::master(comm))
        {
            UIPstream is(UPstream::masterNo(), pBufs);
            is >> data;

            realIsPtr.reset(new IListStream(std::move(data)));
            realIsPtr->name() = fName;
        }
    }

    Pstream::scatter(ok, Pstream::msgType(), comm);

    //- Set stream properties from realIsPtr on master

    // Scatter master header info
    int verValue;
    int fmtValue;
    unsigned streamLabelByteSize;
    unsigned streamScalarByteSize;
    if (UPstream::master(comm))
    {
        verValue = realIsPtr().version().canonical();
        fmtValue = static_cast<int>(realIsPtr().format());
        streamLabelByteSize = realIsPtr().labelByteSize();
        streamScalarByteSize = realIsPtr().scalarByteSize();
    }
    Pstream::scatter(verValue); //,  Pstream::msgType(), comm);
    Pstream::scatter(fmtValue); //,  Pstream::msgType(), comm);
    Pstream::scatter(streamLabelByteSize); //,  Pstream::msgType(), comm);
    Pstream::scatter(streamScalarByteSize); //,  Pstream::msgType(), comm);

    realIsPtr().version(IOstreamOption::versionNumber::canonical(verValue));
    realIsPtr().format(IOstreamOption::streamFormat(fmtValue));
    realIsPtr().setLabelByteSize(streamLabelByteSize);
    realIsPtr().setScalarByteSize(streamScalarByteSize);

    word name(headerIO.name());
    Pstream::scatter(name, Pstream::msgType(), comm);
    headerIO.rename(name);
    Pstream::scatter(headerIO.headerClassName(), Pstream::msgType(), comm);
    Pstream::scatter(headerIO.note(), Pstream::msgType(), comm);
    //Pstream::scatter(headerIO.instance(), Pstream::msgType(), comm);
    //Pstream::scatter(headerIO.local(), Pstream::msgType(), comm);

    return realIsPtr;
}


void Foam::decomposedBlockData::gather
(
    const label comm,
    const label data,
    labelList& datas
)
{
    const label nProcs = UPstream::nProcs(comm);
    datas.setSize(nProcs);

    char* data0Ptr = reinterpret_cast<char*>(datas.data());

    List<int> recvOffsets;
    List<int> recvSizes;
    if (UPstream::master(comm))
    {
        recvOffsets.setSize(nProcs);
        forAll(recvOffsets, proci)
        {
            // Note: truncating long int to int since UPstream::gather limited
            // to ints
            recvOffsets[proci] =
                int(reinterpret_cast<char*>(&datas[proci]) - data0Ptr);
        }
        recvSizes.setSize(nProcs, sizeof(label));
    }

    UPstream::gather
    (
        reinterpret_cast<const char*>(&data),
        sizeof(label),
        data0Ptr,
        recvSizes,
        recvOffsets,
        comm
    );
}


void Foam::decomposedBlockData::gatherSlaveData
(
    const label comm,
    const UList<char>& data,
    const labelUList& recvSizes,

    const label startProc,
    const label nProcs,

    List<int>& sliceOffsets,
    List<char>& recvData
)
{
    // Calculate master data
    List<int> sliceSizes;
    if (UPstream::master(comm))
    {
        const label numProcs = UPstream::nProcs(comm);

        sliceSizes.resize(numProcs, 0);
        sliceOffsets.resize(numProcs+1, 0);

        int totalSize = 0;
        label proci = startProc;
        for (label i = 0; i < nProcs; i++)
        {
            sliceSizes[proci] = int(recvSizes[proci]);
            sliceOffsets[proci] = totalSize;
            totalSize += sliceSizes[proci];
            ++proci;
        }
        sliceOffsets[proci] = totalSize;
        recvData.setSize(totalSize);
    }

    int nSend = 0;
    if
    (
       !UPstream::master(comm)
     && (UPstream::myProcNo(comm) >= startProc)
     && (UPstream::myProcNo(comm) < startProc+nProcs)
    )
    {
        // Note: UPstream::gather limited to int
        nSend = int(data.size_bytes());
    }

    UPstream::gather
    (
        data.cdata(),
        nSend,

        recvData.data(),
        sliceSizes,
        sliceOffsets,
        comm
    );
}


Foam::label Foam::decomposedBlockData::calcNumProcs
(
    const label comm,
    const off_t maxBufferSize,
    const labelUList& recvSizes,
    const label startProci
)
{
    const label nProcs = UPstream::nProcs(comm);

    label nSendProcs = -1;
    if (UPstream::master(comm))
    {
        off_t totalSize = recvSizes[startProci];
        label proci = startProci+1;
        while (proci < nProcs && (totalSize+recvSizes[proci] < maxBufferSize))
        {
            totalSize += recvSizes[proci];
            proci++;
        }

        nSendProcs = proci-startProci;
    }

    // Scatter nSendProcs
    label n;
    UPstream::scatter
    (
        reinterpret_cast<const char*>(&nSendProcs),
        List<int>(nProcs, sizeof(nSendProcs)),
        List<int>(nProcs, Zero),
        reinterpret_cast<char*>(&n),
        sizeof(n),
        comm
    );

    return n;
}


bool Foam::decomposedBlockData::writeBlocks
(
    const label comm,
    autoPtr<OSstream>& osPtr,
    List<std::streamoff>& blockOffset,
    const UList<char>& masterData,

    const labelUList& recvSizes,
    const PtrList<SubList<char>>& slaveData,

    const UPstream::commsTypes commsType,
    const bool syncReturnState,
    const bool separateHeader
)
{
    if (debug)
    {
        Pout<< "decomposedBlockData::writeBlocks:"
            << " stream:" << (osPtr ? osPtr->name() : "invalid")
            << " data:" << masterData.size()
            << " (master only) slaveData:" << slaveData.size()
            << " commsType:" << Pstream::commsTypeNames[commsType] << endl;
    }

    const label nProcs = UPstream::nProcs(comm);

    bool ok = true;

    // Write master data
    if (UPstream::master(comm))
    {
        blockOffset.resize(nProcs);

        OSstream& os = *osPtr;

        label leadingChars = 0;

        if (separateHeader)
        {
            // Check for header in first block
            UIListStream headerStream(masterData);
            token tok(headerStream);

            if (tok.isWord("FoamFile"))
            {
                const char* blockBegin =
                    (masterData.cdata() + headerStream.pos());

                // Read (and discard) sub-dictionary content
                dictionary headerDict(headerStream, false);

                const char* blockEnd =
                    (masterData.cdata() + headerStream.pos());

                leadingChars = label(headerStream.pos());
                while (std::isspace(masterData[leadingChars]))
                {
                    ++leadingChars;
                }

                // Write the header content as FoamBlock
                os.write(word(decomposedBlockData::headerName));

                if (isA<OFstream>(os))
                {
                    // Serial file output - can use writeRaw()
                    os.writeRaw(blockBegin, (blockEnd-blockBegin));
                }
                else
                {
                    // Other cases are less fortunate, and no std::string_view
                    std::string str(blockBegin, (blockEnd-blockBegin));
                    os.writeQuoted(str, false);
                }
                os << nl;
            }
        }

        blockOffset[UPstream::masterNo()] =
            decomposedBlockData::writeBlockEntry
            (
                os,
                UPstream::masterNo(),
                SubList<char>
                (
                    masterData,
                    masterData.size()-leadingChars,
                    leadingChars
                )
            );

        ok = os.good();
    }

    if (slaveData.size())
    {
        // Already have gathered the slave data.

        if (UPstream::master(comm))
        {
            // Master data already written ...
            OSstream& os = *osPtr;

            // Write slaves
            for (label proci = 1; proci < nProcs; ++proci)
            {
                blockOffset[proci] =
                    decomposedBlockData::writeBlockEntry
                    (
                        os,
                        proci,
                        slaveData[proci]
                    );
            }

            ok = os.good();
        }
    }
    else if (commsType == UPstream::commsTypes::scheduled)
    {
        if (UPstream::master(comm))
        {
            // Master data already written ...
            OSstream& os = *osPtr;

            // Receive and write slaves
            DynamicList<char> elems;
            for (label proci = 1; proci < nProcs; ++proci)
            {
                elems.resize(recvSizes[proci]);
                IPstream::read
                (
                    UPstream::commsTypes::scheduled,
                    proci,
                    elems.data(),
                    elems.size_bytes(),
                    Pstream::msgType(),
                    comm
                );

                blockOffset[proci] =
                    decomposedBlockData::writeBlockEntry
                    (
                        os,
                        proci,
                        elems
                    );
            }

            ok = os.good();
        }
        else
        {
            UOPstream::write
            (
                UPstream::commsTypes::scheduled,
                UPstream::masterNo(),
                masterData.cdata(),
                masterData.size_bytes(),
                Pstream::msgType(),
                comm
            );
        }
    }
    else
    {
        // Master data already written ...

        // Find out how many processor can be received into
        // maxMasterFileBufferSize

        // Starting slave processor and number of processors
        label startProc = 1;
        label nSendProcs = nProcs-1;

        while (nSendProcs > 0 && startProc < nProcs)
        {
            nSendProcs = calcNumProcs
            (
                comm,
                off_t
                (
                    fileOperations::masterUncollatedFileOperation::
                    maxMasterFileBufferSize
                ),
                recvSizes,
                startProc
            );

            if (nSendProcs == 0)
            {
                break;
            }


            // Gather data from (a slice of) the slaves
            List<int> sliceOffsets;
            List<char> recvData;
            gatherSlaveData
            (
                comm,
                masterData,
                recvSizes,

                startProc,      // startProc,
                nSendProcs,     // nProcs,

                sliceOffsets,
                recvData
            );

            if (UPstream::master(comm))
            {
                OSstream& os = *osPtr;

                // Write slaves
                for
                (
                    label proci = startProc;
                    proci < startProc+nSendProcs;
                    ++proci
                )
                {
                    SubList<char> dataSlice
                    (
                        recvData,
                        sliceOffsets[proci+1]-sliceOffsets[proci],
                        sliceOffsets[proci]
                    );

                    blockOffset[proci] =
                        decomposedBlockData::writeBlockEntry
                        (
                            os,
                            proci,
                            dataSlice
                        );
                }
            }

            startProc += nSendProcs;
        }

        if (UPstream::master(comm))
        {
            ok = osPtr->good();
        }
    }

    if (syncReturnState)
    {
        //- Enable to get synchronised error checking. Is the one that keeps
        //  slaves as slow as the master (which does all the writing)
        Pstream::scatter(ok, Pstream::msgType(), comm);
    }

    return ok;
}


bool Foam::decomposedBlockData::read()
{
    autoPtr<ISstream> isPtr;
    fileName objPath(fileHandler().filePath(false, *this, word::null));

    bool hasBlockHeader = false;
    if (UPstream::master(comm_))
    {
        isPtr.reset(new IFstream(objPath));
        auto& is = *isPtr;
        IOobject::readHeader(is);

        is.fatalCheck(FUNCTION_NAME);
        token tok(is);
        is.fatalCheck(FUNCTION_NAME);

        hasBlockHeader = tok.isWord(decomposedBlockData::headerName);

        if (hasBlockHeader)
        {
            // FoamBlock: header for block data
            dictionary headerDict(is, false);  // Read sub-dictionary content

            IOobject dummyIO(*this);
            streamOption_ = dummyIO.parseHeader(headerDict);
        }
        else
        {
            is.putBack(tok);
        }
    }

    bool ok = readBlocks(comm_, isPtr, contentData_, commsType_);

    // No separate FoamBlock:
    // Read header from first block

    if (!hasBlockHeader && UPstream::master(comm_))
    {
        IOobject dummyIO(*this);

        UIListStream headerStream(contentData_);
        if (dummyIO.readHeader(headerStream))
        {
            streamOption_ = static_cast<IOstreamOption>(headerStream);
        }
        else
        {
            FatalIOErrorInFunction(headerStream)
                << "Problem while reading header for object "
                << isPtr->name() << nl
                << exit(FatalIOError);
        }
    }

    // Scatter IOstreamOption
    int verValue;
    int fmtValue;
    if (UPstream::master(comm_))
    {
        verValue = streamOption_.version().canonical();
        fmtValue = static_cast<int>(streamOption_.format());
    }
    Pstream::scatter(verValue); //,  Pstream::msgType(), comm);
    Pstream::scatter(fmtValue); //,  Pstream::msgType(), comm);

    streamOption_.version(IOstreamOption::versionNumber::canonical(verValue));
    streamOption_.format(IOstreamOption::streamFormat(fmtValue));

    return ok;
}


bool Foam::decomposedBlockData::writeData(Ostream& os) const
{
    IOobject io(*this);
    IOstreamOption streamOpt(streamOption_);

    // Ensure consistent stream format information
    int verValue;
    int fmtValue;

    bool embeddedHeader = false;
    if (UPstream::master(comm_))
    {
        // Check for header in first block
        UIListStream headerStream(contentData_);
        token tok(headerStream);

        embeddedHeader = tok.isWord("FoamFile");
        if (embeddedHeader)
        {
            headerStream.rewind();

            verValue = headerStream.version().canonical();
            fmtValue = static_cast<int>(headerStream.format());
        }
        else
        {
            verValue = streamOption_.version().canonical();
            fmtValue = static_cast<int>(streamOption_.format());
        }
    }

    // Scatter header information
    Pstream::scatter(verValue, Pstream::msgType(), comm_);
    Pstream::scatter(fmtValue, Pstream::msgType(), comm_);

    streamOpt.version(IOstreamOption::versionNumber::canonical(verValue));
    streamOpt.format(IOstreamOption::streamFormat(fmtValue));

    //word masterName(name());
    //Pstream::scatter(masterName, Pstream::msgType(), comm_);

    Pstream::scatter(io.headerClassName(), Pstream::msgType(), comm_);
    Pstream::scatter(io.note(), Pstream::msgType(), comm_);
    //Pstream::scatter(io.instance(), Pstream::msgType(), comm);
    //Pstream::scatter(io.local(), Pstream::msgType(), comm);

    fileName masterLocation(instance()/db().dbDir()/local());
    Pstream::scatter(masterLocation, Pstream::msgType(), comm_);

    if (!embeddedHeader)
    {
        decomposedBlockData::writeHeader
        (
            os,
            streamOpt,
            io.headerClassName(),
            io.note(),
            masterLocation,
            name()
        );
    }

    // Write the character data
    if (isA<OFstream>(os))
    {
        // Serial file output - can use writeRaw()
        os.writeRaw(contentData_.cdata(), contentData_.size_bytes());
    }
    else
    {
        // Other cases are less fortunate, and no std::string_view
        std::string str(contentData_.cdata(), contentData_.size_bytes());
        os.writeQuoted(str, false);
    }

    if (!Pstream::master(comm_))
    {
        IOobject::writeEndDivider(os);
    }

    return os.good();
}


bool Foam::decomposedBlockData::writeObject
(
    IOstreamOption streamOpt,
    const bool valid
) const
{
    // Remember the output format
    const_cast<IOstreamOption&>(streamOption_) = streamOpt;

    autoPtr<OSstream> osPtr;
    if (UPstream::master(comm_))
    {
        // Note: always write binary. These are strings so readable anyway.
        //       They have already be tokenised on the sending side.

        osPtr.reset(new OFstream(objectPath(), IOstreamOption::BINARY));

        // Update meta-data for current state
        static_cast<regIOobject&>
        (
            const_cast<decomposedBlockData&>(*this)
        ).updateMetaData();

        decomposedBlockData::writeHeader
        (
            *osPtr,
            streamOpt,  // streamOpt for data
            static_cast<const IOobject&>(*this)
        );
    }

    labelList recvSizes;
    gather(comm_, label(contentData_.size_bytes()), recvSizes);

    List<std::streamoff> blockOffsets;
    PtrList<SubList<char>> slaveData;  // dummy slave data
    return writeBlocks
    (
        comm_,
        osPtr,
        blockOffsets,
        contentData_,
        recvSizes,
        slaveData,
        commsType_,
        true,   // syncReturnState
        false   // no separateHeader - already done (if required) above
    );
}


// ************************************************************************* //
