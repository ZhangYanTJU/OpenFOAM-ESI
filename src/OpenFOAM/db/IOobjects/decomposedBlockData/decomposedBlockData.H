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
    A secondary \b FoamBlock header specifies the header information
    for the blocks. Each data block is written as a dictionary entry.
    .
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
}
FoamBlock
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       vectorField;
    location    "constant/polyMesh";
    object      points;
}

processor0 List<char>
NCHARS
(...content...);

processor1 List<char>
NCHARS
(...content...);

...
\endverbatim

    The original format (OpenFOAM-v2012 and earlier) is slightly different.
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

#ifndef decomposedBlockData_H
#define decomposedBlockData_H

#include "IOList.H"
#include "regIOobject.H"
#include "UPstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

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

        //- Name for secondary file header
        static constexpr const char* const headerName = "FoamBlock";

        //- Type to use for gather
        const UPstream::commsTypes commsType_;

        //- Communicator for all parallel comms
        const label comm_;

        //- The format/version for the block(s).
        IOstreamOption streamOption_;

        //- The block content
        List<char> contentData_;


    // Protected Member Functions

        //- Helper: determine number of processors whose recvSizes fits
        //- into maxBufferSize
        static label calcNumProcs
        (
            const label comm,
            const off_t maxBufferSize,
            const labelUList& recvSizes,
            const label startProci
        );

        //- Read data into *this. ISstream is only valid on master.
        static bool readBlocks
        (
            const label comm,
            autoPtr<ISstream>& isPtr,
            List<char>& contentChars,
            const UPstream::commsTypes commsType
        );


public:

    // Static Members

        //- The collated block output format (0: old)
        static int collatedBlockFormat_;

        //- Using the old block output format (2012 and older)
        static bool oldBlockFormat() noexcept
        {
            return collatedBlockFormat_ == 0;
        }


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
        //  The master data may or may not contain an embedded FoamFile header
        virtual bool writeData(Ostream& os) const;

        //- Write using stream options
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool valid
        ) const;


    // Helpers

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
            const word& objectName
        );

        //- Helper: write FoamFile IOobject header
        static void writeHeader
        (
            Ostream& os,
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
            const UList<char>& charData
        );

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
        //- data in stream. Note: isPtr is only valid on master.
        static autoPtr<ISstream> readBlocks
        (
            const label comm,
            const fileName& fName,
            autoPtr<ISstream>& isPtr,
            IOobject& headerIO,
            const UPstream::commsTypes commsType
        );

        //- Helper: gather single label. Note: using native Pstream.
        //  datas sized with num procs but undefined contents on
        //  slaves
        static void gather
        (
            const label comm,
            const label data,
            labelList& datas
        );

        //- Helper: gather data from (subset of) slaves.
        //
        //  Returns:
        //  - recvData : received data
        //  - recvOffsets : offset in data. recvOffsets is nProcs+1
        static void gatherSlaveData
        (
            const label comm,
            const UList<char>& data,
            const labelUList& recvSizes,

            const label startProc,
            const label nProcs,

            List<int>& recvOffsets,
            List<char>& recvData
        );

        //- Write *this. Ostream only valid on master.
        //  Returns offsets of processor blocks in blockOffset
        //
        //  \note the separateHeader parameter will attempt to rewrite
        //  the embedded FoamFile information as a separate FoamBlock header
        static bool writeBlocks
        (
            const label comm,
            autoPtr<OSstream>& osPtr,
            List<std::streamoff>& blockOffset,

            const UList<char>& masterData,

            const labelUList& recvSizes,

            // optional slave data (on master)
            const PtrList<SubList<char>>& slaveData,

            const UPstream::commsTypes,
            const bool syncReturnState = true,
            const bool separateHeader = !decomposedBlockData::oldBlockFormat()
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
