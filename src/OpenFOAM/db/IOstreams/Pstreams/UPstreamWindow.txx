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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::UList<Type> Foam::UPstream::Window::allocate
(
    std::streamsize count,
    UPstream::Communicator communicator,
    const bool shared
)
{
    if constexpr (!is_contiguous_v<Type>)
    {
        FatalErrorInFunction
            << "Only contiguous data can be supported!"
            << Foam::abort(FatalError);
        return UList<Type>();
    }
    else
    {
        auto [ptr, len] =
            mpi_win_allocate(count, sizeof(Type), communicator, shared);
        return UList<Type>(reinterpret_cast<Type*>(ptr), len);
    }
}


template<class Type>
Foam::UList<Type> Foam::UPstream::Window::allocate
(
    std::streamsize count,
    const int communicator,
    const bool shared
)
{
    if constexpr (!is_contiguous_v<Type>)
    {
        FatalErrorInFunction
            << "Only contiguous data can be supported!"
            << Foam::abort(FatalError);
        return UList<Type>();
    }
    else
    {
        auto [ptr, len] =
            mpi_win_allocate(count, sizeof(Type), communicator, shared);
        return UList<Type>(reinterpret_cast<Type*>(ptr), len);
    }
}


template<class Type>
Foam::UList<Type> Foam::UPstream::Window::allocate_shared
(
    std::streamsize count,
    UPstream::Communicator communicator
)
{
    return allocate<Type>(count, communicator, true);
}


template<class Type>
Foam::UList<Type> Foam::UPstream::Window::allocate_shared
(
    std::streamsize count,
    const int communicator
)
{
    return allocate<Type>(count, communicator, true);
}


template<class Type>
bool Foam::UPstream::Window::create
(
    const Type* buffer,
    std::streamsize count,
    UPstream::Communicator communicator
)
{
    if constexpr (!is_contiguous_v<Type>)
    {
        FatalErrorInFunction
            << "Only contiguous data can be supported!"
            << Foam::abort(FatalError);
        return false;
    }
    else
    {
        // const_cast since we cannot specify readonly or read/write
        return mpi_win_create
        (
            const_cast<Type*>(buffer),
            count,
            sizeof(Type),
            communicator
        );
    }
}


template<class Type>
bool Foam::UPstream::Window::create
(
    const Type* buffer,
    std::streamsize count,
    const int communicator
)
{
    if constexpr (!is_contiguous_v<Type>)
    {
        FatalErrorInFunction
            << "Only contiguous data can be supported!"
            << Foam::abort(FatalError);
        return false;
    }
    else
    {
        // const_cast since we cannot specify readonly or read/write
        return mpi_win_create
        (
            const_cast<Type*>(buffer),
            count,
            sizeof(Type),
            communicator
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::UList<Type> Foam::UPstream::Window::view() const
{
    if constexpr (!is_contiguous_v<Type>)
    {
        FatalErrorInFunction
            << "Only contiguous data can be supported!"
            << Foam::abort(FatalError);
        return UList<Type>();
    }
    else
    {
        auto [ptr, len] = mpi_win_query(*this, sizeof(Type));
        return UList<Type>(reinterpret_cast<Type*>(ptr), len);
    }
}


template<class Type>
Foam::UList<Type> Foam::UPstream::Window::view_shared(int target_rank) const
{
    if constexpr (!is_contiguous_v<Type>)
    {
        FatalErrorInFunction
            << "Only contiguous data can be supported!"
            << Foam::abort(FatalError);
        return UList<Type>();
    }
    else
    {
        auto [ptr, count] =
            mpi_win_query_shared(*this, target_rank, sizeof(Type));
        return UList<Type>(reinterpret_cast<Type*>(ptr), count);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::UPstream::Window::get
(
    Type* buffer,
    std::streamsize count,
    int fromProcNo,
    int target_disp
) const
{
    if constexpr (!is_contiguous_v<Type>)
    {
        FatalErrorInFunction
            << "Only contiguous data can be supported!"
            << Foam::abort(FatalError);
        return false;
    }
    else
    {
        // Use element or component type (or byte-wise) for data type
        return this->get_data
        (
            buffer,         // The data or cmpt pointer
            UPstream_dataType<Type>::size(count),
            UPstream_dataType<Type>::datatype_id,
            fromProcNo,
            target_disp
        );
    }
}


template<class Type>
bool Foam::UPstream::Window::put
(
    const Type* buffer,
    std::streamsize count,
    int toProcNo,
    int target_disp
) const
{
    if constexpr (!is_contiguous_v<Type>)
    {
        FatalErrorInFunction
            << "Only contiguous data can be supported!"
            << Foam::abort(FatalError);
        return false;
    }
    else
    {
        // Use element or component type (or byte-wise) for data type
        return this-put_data
        (
            buffer,  // The data or cmpt pointer
            UPstream_dataType<Type>::size(count),
            UPstream_dataType<Type>::datatype_id,
            toProcNo,
            target_disp
        );
    }
}


template<class Type>
bool Foam::UPstream::Window::put
(
    const UPstream::opCodes opCodeId,
    const Type* buffer,
    std::streamsize count,
    int toProcNo,
    int target_disp
) const
{
    if constexpr (!UPstream_basic_dataType<Type>::value)
    {
        FatalErrorInFunction
            << "Only basic data types are supported!"
            << Foam::abort(FatalError);
        return false;
    }
    else
    {
        // Basic types (no user types) only
        return this->put_data
        (
            opCodeId,
            buffer,  // The data or cmpt pointer
            UPstream_basic_dataType<Type>::size(count),
            UPstream_basic_dataType<Type>::datatype_id,
            toProcNo,
            target_disp
        );
    }
}


template<class Type>
bool Foam::UPstream::Window::get
(
    UList<Type>& buffer,
    int fromProcNo,
    int target_disp
) const
{
    return this->get
    (
        buffer.data(), buffer.size(),
        fromProcNo,
        target_disp
    );
}


template<class Type>
bool Foam::UPstream::Window::put
(
    const UList<Type>& buffer,
    int toProcNo,
    int target_disp
) const
{
    return this->put
    (
        buffer.cdata(), buffer.size(),
        toProcNo,
        target_disp
    );
}


template<class Type>
bool Foam::UPstream::Window::put
(
    const UPstream::opCodes opCodeId,
    const UList<Type>& buffer,
    int toProcNo,
    int target_disp
) const
{
    return this->put
    (
        opCodeId,
        buffer.cdata(), buffer.size(),
        toProcNo,
        target_disp
    );
}


template<class Type>
bool Foam::UPstream::Window::get
(
    SubList<Type> buffer,
    int fromProcNo,
    int target_disp
) const
{
    return this->get
    (
         buffer.data(), buffer.size(),
         fromProcNo,
         target_disp
     );
}


template<class Type>
bool Foam::UPstream::Window::fetch_and_op
(
    const UPstream::opCodes opCodeId,
    const Type& origin,
    Type& result,
    int target_rank,
    int target_disp
) const
{
    if constexpr (!UPstream_basic_dataType<Type>::value)
    {
        FatalErrorInFunction
            << "Only basic data types are supported!"
            << Foam::abort(FatalError);
        return false;
    }
    else
    {
        // Basic types (no user types) and a single element only!
        return this->mpi_fetch_and_op
        (
            opCodeId,
           &origin,
           &result,
            UPstream_basic_dataType<Type>::datatype_id,
            target_rank,
            target_disp
        );
    }
}


template<class Type>
bool Foam::UPstream::Window::put
(
    const SubList<Type> buffer,
    int toProcNo,
    int target_disp
) const
{
    return this->put
    (
        buffer.cdata(), buffer.size(),
        toProcNo,
        target_disp
    );
}


template<class Type>
bool Foam::UPstream::Window::put
(
    const UPstream::opCodes opCodeId,
    const SubList<Type> buffer,
    int toProcNo,
    int target_disp
) const
{
    return this->put
    (
        opCodeId,
        buffer.cdata(), buffer.size(),
        toProcNo,
        target_disp
    );
}


// ************************************************************************* //
