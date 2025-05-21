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

#include "argList.H"
#include "IOstreams.H"

// Enable/disable based on header
#ifdef Foam_MemoryPool_H
#define FOAM_HAS_MEMORY_POOL
#else
#undef FOAM_HAS_MEMORY_POOL
#endif

// options
int min_align_size = 5;
int min_pool_size = 10;
bool use_aligned_alloc(true);
bool use_aligned_dealloc(true);

//- True if size exceeds the min-size for using memory alignment
template<class IntType>
inline bool use_alignment(IntType n) noexcept
{
    return (n >= min_align_size);
}


//- True if size exceeds the min-size for using the memory pool
template<class IntType>
inline bool use_memory_pool(IntType n) noexcept
{
    return (n >= IntType(min_pool_size));
}


//- Default alignment
inline constexpr std::align_val_t default_alignment() noexcept
{
    return std::align_val_t(64);
}


//- Allocate from memory pool (if active) or aligned/normal
template<class T, class IntType>
inline T* my_allocate(IntType n)
{
    std::cerr<< "my_allocate(" << n << ")\n";

    if (use_aligned_alloc && use_alignment(n))
    {
        #ifdef FOAM_HAS_MEMORY_POOL
        if
        (
            void *pool_ptr
            (
                // Consider memory pool for large amounts of data
                use_memory_pool(n)
              ? Foam::MemoryPool::try_allocate(sizeof(T)*n)
              : nullptr
            );
            pool_ptr
        )
        {
            // Placement new
            return new (pool_ptr) T[n];
        }
        else
        #endif
        {
            return new (default_alignment()) T[n];
        }
    }
    else
    {
        // Plain new
        return new T[n];
    }
}


//- Deallocate from memory pool or normal
template<class T>
inline void my_deallocate(T* ptr)
{
    std::cerr<< "my_deallocate() : " << Foam::name(ptr) << '\n';

    #ifdef FOAM_HAS_MEMORY_POOL
    if
    (
        ptr && !Foam::MemoryPool::try_deallocate(ptr)
    )
    #endif
    {
        // Plain new
        delete[] ptr;
    }
}


//- Deallocate from memory pool or aligned/normal
template<class T, class IntType>
inline void my_deallocate(T* ptr, IntType n)
{
    std::cerr<< "my_deallocate(" << n << ") : " << Foam::name(ptr) << '\n';

    if (use_aligned_dealloc && use_alignment(n))
    {
        if
        (
            #ifdef FOAM_HAS_MEMORY_POOL
            ptr && !Foam::MemoryPool::try_deallocate(ptr)
            #else
            ptr
            #endif
        )
        {
            ::operator delete[](ptr, default_alignment());
        }
    }
    else
    {
        // Plain new
        delete[] ptr;
    }
}


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();
    argList::addOption
    (
        "min-align",
        "INT",
        "Min number of elements for memory alignment (default: 5)"
    );
    argList::addOption
    (
        "min-pool",
        "INT",
        "Min number of elements for using memory pool (default: 10)"
    );
    argList::addOption
    (
        "count",
        "INT",
        "Number of elements to test (default: 10)"
    );
    argList::addBoolOption
    (
        "no-align",
        "Disable aligned alloc/dealloc"
    );
    argList::addBoolOption
    (
        "no-align-alloc",
        "Disable aligned alloc (default: false)"
    );
    argList::addBoolOption
    (
        "no-align-dealloc",
        "Disable aligned dealloc (default: false)"
    );

    #include "setRootCase.H"

    label count(10);

    args.readIfPresent("count", count);
    args.readIfPresent("min-align", min_align_size);
    args.readIfPresent("min-pool", min_pool_size);

    if (min_pool_size < min_align_size)
    {
        min_pool_size = min_align_size;
    }

    if (args.found("no-align"))
    {
        use_aligned_alloc = false;
        use_aligned_dealloc = false;
    }
    else
    {
        use_aligned_alloc = !args.found("no-align-alloc");
        use_aligned_dealloc = !args.found("no-align-dealloc");
    }


    Info<< "Testing with " << count << " elements" << nl
        << "min-align: " << int(min_align_size) << " elements" << nl
        #ifdef FOAM_HAS_MEMORY_POOL
        << "min-pool: " << int(min_pool_size)
        << " elements, active:" << MemoryPool::active() << nl
        #endif
        << "alignment: " << int(default_alignment()) << " bytes" << nl
        << nl;

    {
        using T = double;
        label len = count;

        UList<T> list(my_allocate<T>(len), len);

        Info<< "List ptr: " << Foam::name(list.data()) << nl;

        list = 1.234;

        Info<< "List: " << list << nl;

        my_deallocate(list.data(), len);

        list = UList<T>();

        my_deallocate(list.data());
        my_deallocate(list.data(), len);
    }

    return 0;
}


// ************************************************************************* //
