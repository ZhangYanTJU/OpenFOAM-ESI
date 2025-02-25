/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
bool Foam::UOPstream::write
(
    const UPstream::commsTypes commsType,
    const int toProcNo,
    const Type* buffer,
    std::streamsize count,
    const int tag,
    const int communicator,
    [[maybe_unused]] UPstream::Request* req,
    const UPstream::sendModes sendMode
)
{
    if constexpr (!is_contiguous_v<Type>)
    {
        // Report parameters to silence compiler warnings about unused
        FatalErrorInFunction
            << "Invalid for non-contiguous data types. "
            << int(commsType) << ':' << int(sendMode)
            << ':' << toProcNo
            << ':' << (buffer != nullptr)
            << ':' << count
            << ':' << tag
            << ':' << communicator
            << Foam::abort(FatalError);
        return false;
    }
    else
    {
        // Use element or component type (or byte-wise) for data type
        return UPstream::mpi_send
        (
            commsType,
            buffer,         // The data or cmpt pointer
            UPstream_dataType<Type>::size(count),
            UPstream_dataType<Type>::datatype_id,
            toProcNo,
            tag,
            communicator,
            req,
            sendMode
       );
    }
}


template<class Type>
bool Foam::UOPstream::write
(
    UPstream::Request& req,
    const int toProcNo,
    const Type* buffer,
    std::streamsize count,
    const int tag,
    const int communicator,
    const UPstream::sendModes sendMode
)
{
    return UOPstream::write
    (
        UPstream::commsTypes::nonBlocking,
        toProcNo,
        buffer,
        count,
        tag,
        communicator,
       &req,
        sendMode
    );
}


template<class Type>
bool Foam::UOPstream::write
(
    const UPstream::commsTypes commsType,
    const int toProcNo,
    const UList<Type>& buffer,
    const int tag,
    const int communicator,
    UPstream::Request* req,
    const UPstream::sendModes sendMode
)
{
    return UOPstream::write
    (
        commsType,
        toProcNo,
        buffer.cdata(), buffer.size(),
        tag,
        communicator,
        req,
        sendMode
    );
}


template<class Type>
bool Foam::UOPstream::write
(
    const UPstream::commsTypes commsType,
    const int toProcNo,
    const SubList<Type> buffer,
    const int tag,
    const int communicator,
    UPstream::Request* req,
    const UPstream::sendModes sendMode
)
{
    return UOPstream::write
    (
        commsType,
        toProcNo,
        buffer.cdata(), buffer.size(),
        tag,
        communicator,
        req,
        sendMode
    );
}


template<class Type>
bool Foam::UOPstream::write
(
    UPstream::Request& req,
    const int toProcNo,
    const UList<Type>& buffer,
    const int tag,
    const int communicator,
    const UPstream::sendModes sendMode
)
{
    return UOPstream::write
    (
        UPstream::commsTypes::nonBlocking,
        toProcNo,
        buffer.cdata(), buffer.size(),
        tag,
        communicator,
       &req,
        sendMode
    );
}


template<class Type>
bool Foam::UOPstream::write
(
    UPstream::Request& req,
    const int toProcNo,
    const SubList<Type> buffer,
    const int tag,
    const int communicator,
    const UPstream::sendModes sendMode
)
{
    return UOPstream::write
    (
        UPstream::commsTypes::nonBlocking,
        toProcNo,
        buffer.cdata(), buffer.size(),
        tag,
        communicator,
       &req,
        sendMode
    );
}


// ************************************************************************* //
