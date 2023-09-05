/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenFOAM Foundation
    Copyright (C) 2020-2025 OpenCFD Ltd.
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
    Foam::decomposedBlockData

Description
    The decomposedBlockData comprise a \c List\<char\> for each output
    processor, typically with IO on the master processor only.

    For decomposedBlockData, we make a distinction between the container
    description and the individual block contents.

    The \b FoamFile header specifies the container characteristics and thus
    has \c class = \c %decomposedBlockData and normally \c format = \c binary.
    This description refers to the \em entire file container, not the
    individual blocks.

    Each processor block is simply a binary chunk of characters and the
    first block also contains the header description for all of the blocks.
    For example,
\verbatim
FoamFile
{
    version     2.0;
    format      binary;
    arch        "LSB;label=32;scalar=64";
    class       decomposedBlockData;
    location    "constant/polyMesh";
    object      points;
    data.format ascii;          // optional
    data.class  vectorField;    // optional
}

// processor0
NCHARS
(FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       vectorField;
    location    "constant/polyMesh";
    object      points;
}
...content...
)

// processor1
NCHARS
(...content...)

...
\endverbatim


SourceFiles
    decomposedBlockData.C
    decomposedBlockDataHeader.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_decomposedBlockData_H
#define Foam_decomposedBlockData_H

#include "regIOobject.H"
#include "ISstream.H"
#include "OSstream.H"
#include "UPstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;

/*---------------------------------------------------------------------------*\
                     Class decomposedBlockData Declaration
\*---------------------------------------------------------------------------*/

class decomposedBlockData
:
    public regIOobject
{
    // Private Functions

        //- Helper: write content for FoamFile IOobject header
        static void writeHeaderContent
        (
            Ostream& os,
            IOstreamOption streamOptContainer,
            const word& objectType,
            const string& note,
            const fileName& location,
            const word& objectName
        );


protected:

    // Protected Data

        //- Type to use for gather
        const UPstream::commsTypes commsType_;

        //- Communicator for all parallel comms
        const label comm_;

        //- The block content
        List<char> contentData_;


    // Protected Member Functions

        //- Read data (on master) and transmit.
        static bool readBlocks
        (
            const label comm,
            // [in] The input stream (only valid on master)
            autoPtr<ISstream>& isPtr,
            // [out] The processor local data
            List<char>& localData,
            const UPstream::commsTypes commsType  /* unused */
        );

        //- Helper: skip a block of (binary) character data
        static bool skipBlockEntry(Istream& is);

public:

    //- Declare type-name, virtual type (with debug switch)
    TypeName("decomposedBlockData");


    // Constructors

        //- Construct given an IOobject
        decomposedBlockData
        (
            const label comm,
            const IOobject& io,
            const UPstream::commsTypes = UPstream::commsTypes::scheduled
        );


    //- Destructor
    virtual ~decomposedBlockData() = default;


    // Member Functions

        //- Read object
        virtual bool read();

        //- Write separated content (assumes content is the serialised data)
        //  The serialised master data should also contain a FoamFile header
        virtual bool writeData(Ostream& os) const;

        //- Write using stream options
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool writeOnProc
        ) const;


    // Helpers

        //- True if object type is a known collated type
        static bool isCollatedType(const word& objectType);

        //- True if object header class is a known collated type
        static bool isCollatedType(const IOobject& io);

        //- Extract number of decomposedBlockData block entries, optionally
        //- with an upper limit.
        //- The input stream should be in a rewound state
        //- (or only have read the header) before calling.
        static label getNumBlocks(Istream& is, const label maxNumBlocks = -1);

        //- True if the given block number (starts at 0) has a corresponding
        //- decomposedBlockData block entry.
        //- The input stream should be in a rewound state
        //- (or only have read the header) before calling.
        //  This will be faster than checking against getNumBlocks()
        //  since it can potentially exit without scanning the entire file.
        static bool hasBlock(Istream& is, const label blockNumber);

        //- Read header as per IOobject with additional handling of
        //- decomposedBlockData
        static bool readHeader(IOobject& io, Istream& is);

        //- Helper: write FoamFile IOobject header
        static void writeHeader
        (
            Ostream& os,
            IOstreamOption streamOptContainer,
            const word& objectType,
            const string& note,
            const fileName& location,
            const word& objectName,
            const dictionary& extraEntries
        );

        //- Helper: write FoamFile IOobject header
        static void writeHeader
        (
            Ostream& os,
            IOstreamOption streamOptData,
            const IOobject& io
        );

        //- Helper: generate additional entries for FoamFile header
        static void writeExtraHeaderContent
        (
            dictionary& dict,
            IOstreamOption streamOptData,
            const IOobject& io
        );

        //- Helper: read block of (binary) character data
        static bool readBlockEntry
        (
            Istream& is,
            List<char>& charData
        );

        //- Helper: write block of (binary) character data
        static std::streamoff writeBlockEntry
        (
            OSstream& os,
            const label blocki,
            const char* str,
            const size_t len
        );

        //- Helper: write block of (binary) character data
        static std::streamoff writeBlockEntry
        (
            OSstream& os,
            const label blocki,
            const UList<char>& s
        )
        {
            return writeBlockEntry(os, blocki, s.cdata(), s.size_bytes());
        }

        //- Helper: write block of (binary) character content
        //  Housekeeping
        static std::streamoff writeBlockEntry
        (
            OSstream& os,
            const label blocki,
            std::string_view sv
        )
        {
            return writeBlockEntry(os, blocki, sv.data(), sv.size());
        }

        //- Helper: write block of (binary) character data
        //  \return -1 on error
        static std::streamoff writeBlockEntry
        (
            OSstream& os,
            IOstreamOption streamOptData,
            const regIOobject& io,
            const label blocki,
            const bool withLocalHeader
        );

        //- Read selected block (non-seeking) + header information
        static autoPtr<ISstream> readBlock
        (
            const label blocki,
            ISstream& is,
            IOobject& headerIO
        );

        //- Read master header information (into headerIO) and return
        //- data in stream.
        static autoPtr<ISstream> readBlocks
        (
            const label comm,
            const fileName& fName,
            //! [in] The input stream (only valid on master)
            autoPtr<ISstream>& isPtr,
            //! [out] header information
            IOobject& headerIO,
            const UPstream::commsTypes commsType  /* unused */
        );

        //- Helper: gather data from (subset of) sub-ranks.
        //  In non-blocking mode it sets up send/recv for non-empty content.
        //  In blocking/scheduled mode it uses MPI_Gatherv to collect data.
        //
        //  Returns:
        //  - recvData : the received data
        //  - recvOffsets : offset in data. recvOffsets is nProcs+1
        static void gatherProcData
        (
            const label comm,
            const UList<char>& localData,   //!< [in] required on all procs
            const labelUList& recvSizes,    //!< [in] only required on master

            const labelRange& whichProcs,   //!< [in] required on all procs

            List<int>& recvOffsets,         //!< [out] only relevant on master
            DynamicList<char>& recvData,    //!< [out] only relevant on master

            const UPstream::commsTypes commsType
        );

        //- Write *this. Ostream only valid on master.
        //  Returns offsets of processor blocks in blockOffset
        static bool writeBlocks
        (
            const label comm,

            //! [in] output stream (relevant on master)
            autoPtr<OSstream>& osPtr,
            //! [out] start offsets to each block (relevant on master),
            //! ignored if List::null() type
            List<std::streamoff>& blockOffset,

            const UList<char>& localData,   //!< [in] required on all procs
            const labelUList& recvSizes,    //!< [in] only required on master

            //! Optional proc data (only written on master)
            //! but \b must also be symmetrically defined (empty/non-empty)
            //! on all ranks
            const UList<std::string_view>& procData,

            const UPstream::commsTypes commsType,
            const bool syncReturnState = true
        );


    // Housekeeping

        //- Write *this. Ostream only valid on master.
        //  Returns offsets of processor blocks in blockOffset
        FOAM_DEPRECATED_FOR(2023-09, "write with std::string_view instead")
        static bool writeBlocks
        (
            const label comm,
            autoPtr<OSstream>& osPtr,
            List<std::streamoff>& blockOffset,

            const UList<char>& localData, // [in] required on all procs
            const labelUList& recvSizes,  // [in] only required on master

            // Optional proc data (only written on master)
            // but \b must also be symmetrically defined (empty/non-empty)
            // on all ranks
            const UPtrList<SubList<char>>& procData,

            const UPstream::commsTypes commsType,
            const bool syncReturnState = true
        )
        {
            // Transcribe to string_view
            List<std::string_view> spans(procData.size());
            forAll(procData, proci)
            {
                const auto* ptr = procData.get(proci);
                if (ptr && !ptr->empty())
                {
                    spans[proci] = std::string_view(ptr->cdata(), ptr->size());
                }
            }

            bool ok = decomposedBlockData::writeBlocks
            (
                comm,
                osPtr,
                blockOffset,
                localData,
                recvSizes,
                spans,
                commsType,
                syncReturnState
            );
            return ok;
        }

        //- Deprecated(2023-09) - consider UPstream::listGatherValue
        //  The only difference is that this gather also resizes the output
        //  on the non-master procs
        //  \deprecated(2023-09) - consider UPstream::listGatherValue
        FOAM_DEPRECATED_FOR(2023-09, "consider UPstream::listGatherValue()")
        static void gather
        (
            const label comm,
            const label localValue,
            labelList& allValues
        )
        {
            allValues.resize_nocopy(UPstream::nProcs(comm));

            UPstream::mpiGather
            (
                reinterpret_cast<const char*>(&localValue),
                allValues.data_bytes(),
                sizeof(label),  // The send/recv size per rank
                comm
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
