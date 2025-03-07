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

Class
    Foam::UPstream::Window

Description
    An opaque wrapper for MPI_Win with a vendor-independent
    representation and without any \c <mpi.h> header dependency.

    The handling for window declaration is very generous - it does not
    distinguish between readonly and read/write windows. This may become
    more restrictive in the future.

Note
    The MPI standard states that MPI_Win is always an opaque object.
    Generally it is either an integer (eg, mpich) or a pointer (eg, openmpi).

SourceFiles
    UPstreamWindow.txx

\*---------------------------------------------------------------------------*/

#ifndef Foam_UPstreamWindow_H
#define Foam_UPstreamWindow_H

#include "UPstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class UPstream::Window Declaration
\*---------------------------------------------------------------------------*/

class UPstream::Window
{
public:

    // Public Types

        //- Storage for MPI_Win (as integer or pointer)
        typedef std::intptr_t value_type;


private:

    // Private Data

        //- The MPI_Win (as wrapped value)
        value_type value_;


protected:

    // Protected Method Functions

        //- Sizing helper for disp_unit
        template<class Type>
        static constexpr auto element_width() noexcept
        {
            if constexpr (std::is_void_v<Type>) { return 1; }
            else { return sizeof(Type); }
        }

        //- Allocate a local or shared memory window.
        //- Uses MPI_Win_allocate() or MPI_Win_allocate_shared(), respectively.
        //  A no-op if not running in parallel.
        //  This is a \em collective call.
        //  \returns (address, size) tuple
        std::pair<void*, int64_t> mpi_win_allocate
        (
            //! Number of elements in the window (non-negative integer)
            std::streamsize num_elements,
            //! Local unit size for displacements (positive integer)
            //! == sizeof(Type)
            int disp_unit,
            //! Communicator (wrapped)
            UPstream::Communicator communicator,
            //! Create shared memory
            const bool shared = false
        );

        //- Allocate a local or shared memory window.
        //- Uses MPI_Win_allocate() or MPI_Win_allocate_shared(), respectively.
        //  A no-op if not running in parallel.
        //  This is a \em collective call.
        //  \returns (address, size) tuple
        std::pair<void*, int64_t> mpi_win_allocate
        (
            //! Number of elements in the window (non-negative integer)
            std::streamsize num_elements,
            //! Local unit size for displacements (positive integer)
            //! == sizeof(Type)
            int disp_unit,
            //! Communicator (internally indexed)
            int communicator,
            //! Create shared memory
            const bool shared = false
        );

        //- Create window onto existing memory with MPI_Win_create().
        //  A no-op if not running in parallel.
        //  This is a \em collective call.
        bool mpi_win_create
        (
            //! Address of window
            void *baseptr,
            //! Number of elements in the window (non-negative integer)
            std::streamsize num_elements,
            //! Local unit size for displacements (positive integer)
            //! == sizeof(Type)
            int disp_unit,
            //! Communicator (wrapped)
            UPstream::Communicator communicator
        );

        //- Create window onto existing memory with MPI_Win_create().
        //  A no-op if not running in parallel.
        //  This is a \em collective call.
        bool mpi_win_create
        (
            //! Address of window
            void *baseptr,
            //! Number of elements in the window (non-negative integer)
            std::streamsize num_elements,
            //! Local unit size for displacements (positive integer)
            //! == sizeof(Type)
            int disp_unit,
            //! Communicator (internally indexed)
            int communicator
        );

        //- Retrieve window sizing information as address/count tuple.
        //- The expected sizeof(Type) is supplied as an assertion parameter.
        //
        //  \returns (nullptr, 0) if not running in parallel or if either
        //  the address or count are zero
        static std::pair<void*, int64_t> mpi_win_query
        (
            UPstream::Window window,
            //! Expected element size (fatal for mismatch with queried value)
            const int expected_disp_unit
        );

        //- Retrieve shared window information as address/count tuple.
        //- The expected sizeof(Type) is supplied as an assertion parameter.
        //
        //  \returns (nullptr, 0) if not running in parallel or if either
        //  the address or count are zero
        static std::pair<void*, int64_t> mpi_win_query_shared
        (
            UPstream::Window window,
            //! The rank to query
            int target_rank,
            //! Expected element size (fatal for mismatch with queried value)
            const int expected_disp_unit
        );


    // Protected Member Functions

        //- Get buffer contents from given rank.
        //  A no-op if not running in parallel and for empty or null buffer.
        //
        //  \note The method uses a \c void pointer and the required data type
        //  (as per MPI). This means it should almost never be called directly
        //  but always via a compile-time checked caller.
        //  \return True on success
        bool get_data
        (
            //! Destination buffer
            void* origin,
            //! The data count - number of elements. Identical on both sides.
            std::streamsize count,
            const UPstream::dataTypes dataTypeId,
            int target_rank,
            int target_disp = 0     //!< Window displacement on target
        ) const;

        //- Put buffer contents to given rank.
        //  A no-op if not running in parallel and for empty or null buffer.
        bool put_data
        (
            //! Source buffer
            const void* origin,
            //! The data count - number of elements. Identical on both sides.
            std::streamsize count,
            const UPstream::dataTypes dataTypeId,
            int target_rank,
            int target_disp = 0     //!< Window displacement on target
        ) const;

        //- Put accumulate buffer contents to given rank.
        //  A no-op if not running in parallel and for empty or null buffer.
        bool put_data
        (
            //! The op-code for accumulate
            const UPstream::opCodes opCodeId,
            //! Source buffer
            const void* origin,
            //! The data count - number of elements. Identical on both sides.
            std::streamsize count,
            const UPstream::dataTypes dataTypeId,
            int target_rank,
            int target_disp = 0     //!< Window displacement on target
        ) const;

        //- Retrieve the remote content (a single value) and then combine
        //- in new content.
        //  A no-op if not running in parallel and for empty or null buffer.
        bool mpi_fetch_and_op
        (
            //! The op-code for accumulate
            const UPstream::opCodes opCodeId,
            //! [in] the content to combine with target
            const void* origin,
            //! [out] target content before the operation
            void* result,
            const UPstream::dataTypes dataTypeId,
            int target_rank,
            int target_disp = 0     //!< Window displacement on target
        ) const;

        //- Entry point to MPI_Win_flush(), MPI_Win_flush_all(),
        //- MPI_Win_flush_local(), MPI_Win_flush_local_all().
        //  Uses rank == -1 to signal 'all'
        void mpi_win_flushing(int rank, bool local=false);

        //- Entry point to MPI_Win_lock(), MPI_Win_lock_all(),
        //- optionally as exclusive lock.
        //  Uses rank == -1 to signal 'all'
        void mpi_win_locking(int rank, bool exclusive=false);

        //- Entry point to MPI_Win_unlock(), MPI_Win_unlock_all().
        //  Uses rank == -1 to signal 'all'
        void mpi_win_unlocking(int rank);


public:

    // Generated Methods

        //- Copy construct
        Window(const Window&) noexcept = default;

        //- Move construct
        Window(Window&&) noexcept = default;

        //- Copy assignment
        Window& operator=(const Window&) noexcept = default;

        //- Move assignment
        Window& operator=(Window&&) noexcept = default;


    // Member Operators

        //- Test for equality
        bool operator==(const Window& rhs) const noexcept
        {
            return (value_ == rhs.value_);
        }

        //- Test for inequality
        bool operator!=(const Window& rhs) const noexcept
        {
            return (value_ != rhs.value_);
        }


    // Constructors

        //- Default construct as MPI_WIN_NULL
        Window() noexcept;

        //- Construct from MPI_Win (as pointer type)
        explicit Window(const void* p) noexcept
        :
            value_(reinterpret_cast<value_type>(p))
        {}

        //- Construct from MPI_Win (as integer type)
        explicit Window(value_type val) noexcept
        :
            value_(val)
        {}


    // Member Functions

    // Basic handling

        //- Return raw value
        value_type value() const noexcept { return value_; }

        //- Return as pointer value
        const void* pointer() const noexcept
        {
            return reinterpret_cast<const void*>(value_);
        }

        //- True if not equal to MPI_WIN_NULL
        bool good() const noexcept;

        //- Reset to default constructed value (MPI_WIN_NULL)
        void reset() noexcept;

        //- The number of ranks associated with the window group.
        //  The same as querying the original communicator, assuming the
        //  communicator is available within the current code scope.
        int size() const;


    // Window creation with allocation

        //- Allocate a local memory region and create window onto it.
        //  A no-op if not running in parallel.
        //  This is a \em collective call.
        template<class Type>
        inline UList<Type> allocate
        (
            //! Number of \b elements to allocate
            std::streamsize count,
            UPstream::Communicator communicator,
            //! Create shared memory region
            const bool shared = false
        );

        //- Allocate a local memory region and create window onto it.
        //  A no-op if not running in parallel.
        //  This is a \em collective call.
        template<class Type>
        inline UList<Type> allocate
        (
            //! Number of \b elements to allocate
            std::streamsize count,
            const int communicator,
            //! Create shared memory region
            const bool shared = false
        );

        //- Allocate a shared memory region and create window onto it.
        //  A no-op if not running in parallel.
        //  This is a \em collective call.
        template<class Type>
        inline UList<Type> allocate_shared
        (
            //! Number of \b elements to allocate
            std::streamsize count,
            UPstream::Communicator communicator
        );

        //- Create a window by allocating a new shared memory region.
        //  A no-op if not running in parallel.
        //  This is a \em collective call.
        template<class Type>
        inline UList<Type> allocate_shared
        (
            //! Number of \b elements to allocate
            std::streamsize count,
            const int communicator
        );


    // Window creation with existing memory

        //- A window exposing a zero-sized memory region.
        //  A no-op if not running in parallel.
        //  This is a \em collective call.
        template<class Type = void>
        bool create(std::nullptr_t, UPstream::Communicator comm)
        {
            return mpi_win_create(nullptr, 0, element_width<Type>(), comm);
        }

        //- A window exposing a zero-sized memory region.
        //  A no-op if not running in parallel.
        //  This is a \em collective call.
        template<class Type = void>
        bool create(std::nullptr_t, const int comm)
        {
            return mpi_win_create(nullptr, 0, element_width<Type>(), comm);
        }

        //- A window exposing an existing memory region.
        //  A no-op if not running in parallel.
        //  This is a \em collective call.
        //  \note Ignores constness since we can't specify a priori
        //      if this should be a read/write or read only buffer.
        template<class Type>
        inline bool create
        (
            const Type* buffer,
            //! Number of elements in the window
            std::streamsize count,
            UPstream::Communicator communicator
        );

        //- A window exposing an existing memory region.
        //  A no-op if not running in parallel.
        //  This is a \em collective call.
        template<class Type>
        inline bool create
        (
            const Type* buffer,
            //! Number of elements in the window
            std::streamsize count,
            const int communicator
        );

        //- A window exposing the specified buffer contents.
        //  A no-op if not running in parallel.
        //  This is a \em collective call.
        template<class Type>
        bool create(const UList<Type>& buffer, UPstream::Communicator comm)
        {
            return create(buffer.cdata(), buffer.size(), comm);
        }

        //- A window exposing the specified buffer contents.
        //  A no-op if not running in parallel.
        //  This is a \em collective call.
        template<class Type>
        bool create(const UList<Type>& buffer, int communicator)
        {
            return create(buffer.cdata(), buffer.size(), communicator);
        }

        //- A window exposing the specified buffer contents.
        //  A no-op if not running in parallel.
        //  This is a \em collective call.
        template<class Type>
        bool create(SubList<Type> buffer, UPstream::Communicator communicator)
        {
            return create(buffer.cdata(), buffer.size(), communicator);
        }

        //- A window exposing the specified buffer contents.
        //  A no-op if not running in parallel.
        //  This is a \em collective call.
        template<class Type>
        bool create(SubList<Type> buffer, int communicator)
        {
            return create(buffer.cdata(), buffer.size(), communicator);
        }


    // Synchronization and resource management

        //- MPI_Win_lock() for given target rank (no assertions),
        //- optionally as exclusive lock.
        void lock(int rank, bool exclusive=false)
        {
            mpi_win_locking(rank, exclusive);
        }

        //- MPI_Win_unlock() for given target rank.
        void unlock(int rank) { mpi_win_unlocking(rank); }

        //- MPI_Win_lock_all(), optionally as exclusive lock.
        void lock_all(bool exclusive=false) { mpi_win_locking(-1, exclusive); }

        //- MPI_Win_unlock_all()
        void unlock_all() { mpi_win_unlocking(-1); }

        //- MPI_Win_flush() for given target rank.
        //  Complete all outstanding RMA operations initiated by the
        //  calling process to the target rank.
        void flush(int rank) { mpi_win_flushing(rank); }

        //- MPI_Win_flush_all()
        //  Complete all outstanding RMA operations at both the origin
        //  and the target
        void flush_all() { mpi_win_flushing(-1); }

        //- MPI_Win_flush_local()
        //  Locally complete at the origin all outstanding RMA
        //  operations initiated by the calling process to the target
        //  process specified by rank.
        void flush_local(int rank) { mpi_win_flushing(rank, true); }

        //- MPI_Win_flush_local_all()
        //  Locally complete at the origin all outstanding RMA
        //  operations to all targets.
        void flush_local_all() { mpi_win_flushing(-1, true); }

        //- MPI_Win_sync() - ignored if the window is not active.
        void sync();

        //- MPI_Win_free().
        //- Closes the window view and frees any associated memory,
        //  eg, from allocate() or allocate_shared().
        //  Ignored if the window is not active.
        //  This is a \em collective call.
        void close();


    // Queries

        //- Test if the window is a shared memory window
        bool is_shared(const bool failNonShared = false) const;

        //- Return view of the currently exposed window content.
        //- No restriction on the type of memory associated with the window.
        //  A no-op (empty list) if not running in parallel
        //  or the window is not active.
        template<class Type>
        inline UList<Type> view() const;

        //- Return view of shared memory window content.
        //  A no-op (empty list) if not running in parallel.
        //  Undefined behaviour (likely Fatal) if a shared memory window
        //  has not been allocated.
        template<class Type>
        inline UList<Type> view_shared(int target_rank) const;


    // Window Access

        //- Get buffer contents from given rank.
        //- A no-op for an empty or null buffer, or if not running in parallel.
        template<class Type>
        inline bool get
        (
            //! Destination buffer
            Type* buffer,
            //! The data count - number of elements. Identical on both sides.
            std::streamsize count,
            int fromProcNo,
            int target_disp = 0     //!< Window displacement on target
        ) const;

        //- Put buffer contents to given rank.
        //  A no-op for an empty/null buffer, or if not running in parallel.
        template<class Type>
        inline bool put
        (
            //! Source buffer
            const Type* buffer,
            //! The data count - number of elements. Identical on both sides.
            std::streamsize count,
            int toProcNo,
            int target_disp = 0     //!< Window displacement on target
        ) const;

        //- Put accumulate buffer contents to given rank.
        //  A no-op for an empty/null buffer, or if not running in parallel.
        //  \note Must correspond to basic data types!
        template<class Type>
        inline bool put
        (
            const UPstream::opCodes opCodeId,
            //! Source buffer
            const Type* buffer,
            //! The data count - number of elements. Identical on both sides.
            std::streamsize count,
            int toProcNo,
            int target_disp = 0     //!< Window displacement on target
        ) const;

        //- Get a single value from given rank.
        //  \note Use persistent data (not temporary) for value
        template<class Type>
        inline bool get_value
        (
            Type& value,
            int fromProcNo,
            int target_disp = 0     //!< Window displacement on target
        ) const
        {
            return this->get(&value, 1, fromProcNo, target_disp);
        }

        //- Put a single value to given rank.
        //  \note Use persistent data (not temporary) for value
        template<class Type>
        inline bool put_value
        (
            const Type& value,
            int toProcNo,
            int target_disp = 0     //!< Window displacement on target
        ) const
        {
            return this->put(&value, 1, toProcNo, target_disp);
        }

        //- Put and accumulate a single value to given rank
        //  \note Use persistent data (not temporary) for value
        template<class Type>
        inline bool put_value
        (
            const UPstream::opCodes opCodeId,
            const Type& value,
            int toProcNo,
            int target_disp = 0     //!< Window displacement on target
        ) const
        {
            return this->put(opCodeId, &value, 1, toProcNo, target_disp);
        }

        //- Get into List storage (contiguous data only)
        //- from window location on given processor.
        //- A no-op for an empty or null buffer, or if not running in parallel.
        //  \note Only valid for contiguous data types.
        template<class Type>
        inline bool get
        (
            UList<Type>& buffer,
            //! Offset number of elements on the toProcNo rank
            int fromProcNo,
            int target_disp = 0     //!< Window displacement on target
        ) const;

        //- Put from List storage (contiguous data only)
        //- to window location on given processor.
        //  \note Only valid for contiguous data types.
        template<class Type>
        inline bool put
        (
            const UList<Type>& buffer,
            int toProcNo,
            int target_disp = 0     //!< Window displacement on target
        ) const;

        //- Put and accumulate from List storage (contiguous data only)
        //- to window location on given processor.
        //- A no-op for an empty or null buffer, or if not running in parallel.
        //  \note Only valid for contiguous data types.
        template<class Type>
        inline bool put
        (
            const UPstream::opCodes opCodeId,
            const UList<Type>& buffer,
            int toProcNo,
            int target_disp = 0     //!< Window displacement on target
        ) const;

        //- Get into SubList storage (contiguous data only)
        //- from window location on given processor.
        //- A no-op for an empty or null buffer, or if not running in parallel.
        //  \note Only valid for contiguous data types.
        template<class Type>
        inline bool get
        (
            SubList<Type> buffer,
            int fromProcNo,
            int target_disp = 0     //!< Window displacement on target
        ) const;

        //- Combine the value of origin into the target and return
        //- the resulting value (MPI_Fetch_and_op).
        //  \note Use persistent data (not temporary) for origin and result,
        //  which must also be disjoint locations.
        //  Only valid for single element (no aggregates).
        //  \return the old target value (before the operation)
        template<class Type>
        inline bool fetch_and_op
        (
            const UPstream::opCodes opCodeId,
            //! the content to combine with target
            const Type& origin,
            //! [out] the old target value before the operation
            Type& result,
            int target_rank,        //!< The target rank
            int target_disp = 0     //!< Window displacement on target
        ) const;


        //- Put from SubList storage (contiguous data only)
        //- to window location on given processor.
        //- A no-op for an empty or null buffer, or if not running in parallel.
        //  \note Only valid for contiguous data types.
        template<class Type>
        inline bool put
        (
            const SubList<Type> buffer,
            int toProcNo,
            int target_disp = 0     //!< Window displacement on target
        ) const;

        //- Put and accumulate from SubList storage (contiguous data only)
        //- to window location on given processor.
        //- A no-op for an empty or null buffer, or if not running in parallel.
        //  \note Only valid for contiguous data types.
        template<class Type>
        inline bool put
        (
            const UPstream::opCodes opCodeId,
            const SubList<Type> buffer,
            int toProcNo,
            int target_disp = 0     //!< Window displacement on target
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UPstreamWindow.txx"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
