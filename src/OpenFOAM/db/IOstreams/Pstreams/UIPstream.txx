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
std::streamsize Foam::UIPstream::read
(
    const UPstream::commsTypes commsType,
    const int fromProcNo,
    Type* buffer,
    std::streamsize count,
    const int tag,
    const int communicator,
    [[maybe_unused]] UPstream::Request* req
)
{
    if constexpr (!is_contiguous_v<Type>)
    {
        // Report parameters to silence compiler warnings about unused
        FatalErrorInFunction
            << "Invalid for non-contiguous data types. "
            << int(commsType) << ':' << fromProcNo
            << ':' << (buffer != nullptr)
            << ':' << count
            << ':' << tag
            << ':' << communicator
            << Foam::abort(FatalError);
        return 0;
    }
    else
    {
        // Use element or component type (or byte-wise) for data type
        return UPstream::mpi_receive
        (
            commsType,
            buffer,         // The data or cmpt pointer
            UPstream_dataType<Type>::size(count),
            UPstream_dataType<Type>::datatype_id,
            fromProcNo,
            tag,
            communicator,
            req
        );
    }
}


template<class Type>
std::streamsize Foam::UIPstream::read
(
    UPstream::Request& req,
    const int fromProcNo,
    Type* buffer,
    std::streamsize count,
    const int tag,
    const int communicator
)
{
    return UIPstream::read
    (
        UPstream::commsTypes::nonBlocking,
        fromProcNo,
        buffer,
        count,
        tag,
        communicator,
       &req
    );
}


template<class Type>
std::streamsize Foam::UIPstream::read
(
    const UPstream::commsTypes commsType,
    const int fromProcNo,
    UList<Type>& buffer,
    const int tag,
    const int communicator,
    UPstream::Request* req
)
{
    return UIPstream::read
    (
        commsType,
        fromProcNo,
        buffer.data(), buffer.size(),
        tag,
        communicator,
        req
    );
}


template<class Type>
std::streamsize Foam::UIPstream::read
(
    const UPstream::commsTypes commsType,
    const int fromProcNo,
    SubList<Type> buffer,
    const int tag,
    const int communicator,
    UPstream::Request* req
)
{
    return UIPstream::read
    (
        commsType,
        fromProcNo,
        buffer.data(), buffer.size(),
        tag,
        communicator,
        req
    );
}


template<class Type>
std::streamsize Foam::UIPstream::read
(
    UPstream::Request& req,
    const int fromProcNo,
    UList<Type>& buffer,
    const int tag,
    const int communicator
)
{
    return UIPstream::read
    (
        UPstream::commsTypes::nonBlocking,
        fromProcNo,
        buffer.data(), buffer.size(),
        tag,
        communicator,
        &req
    );
}


template<class Type>
std::streamsize Foam::UIPstream::read
(
    UPstream::Request& req,
    const int fromProcNo,
    SubList<Type> buffer,
    const int tag,
    const int communicator
)
{
    return UIPstream::read
    (
        UPstream::commsTypes::nonBlocking,
        fromProcNo,
        buffer.data(), buffer.size(),
        tag,
        communicator,
       &req
    );
}


// ************************************************************************* //
