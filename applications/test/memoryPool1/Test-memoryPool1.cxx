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
#include "MemoryPool.H"
#include "IOstreams.H"

// options
int min_size = 5;
bool use_aligned_alloc(true);
bool use_aligned_dealloc(true);
bool use_dynamic_alignment(false);

// Some type of selector
template<class IntType>
inline bool use_pool(IntType n) noexcept
{
    return (n >= IntType(min_size));
}


//- Default alignment that is size-independent
inline std::align_val_t constexpr my_align() noexcept
{
    return std::align_val_t(64);
}


//- Default alignment that is size-dependent (depends on number of elements)
template<class IntType>
inline std::align_val_t constexpr my_align(IntType n) noexcept
{
    return std::align_val_t(n > 32 ? 256 : 16);
}


//- Allocate from memory pool (if active) or aligned/normal
template<class T, class IntType>
inline T* my_allocate(IntType n)
{
    std::cerr<< "my_allocate(" << n << ")\n";

    if
    (
        void *pool_ptr
        (
            // Consider memory pool for large amounts of data
            use_pool(n)
          ? Foam::MemoryPool::try_allocate(sizeof(T)*n)
          : nullptr
        );
        pool_ptr
    )
    {
        // Placement new
        return new (pool_ptr) T[n];
    }
    else if (use_aligned_alloc)
    {
        if (!use_dynamic_alignment)
        {
            // Constant alignment
            return new (my_align()) T[n];
        }
        else
        {
            // Alignment depends on the number of elements
            return new (my_align(n)) T[n];
        }
    }
    else
    {
        // Plain new
        return new T[n];
    }
}


//- Deallocate from memory pool or aligned/normal
template<class T, class IntType>
inline void my_deallocate(T* ptr, IntType n)
{
    std::cerr<< "my_deallocate(" << n << ")\n";

    if (ptr)
    {
        if (!Foam::MemoryPool::try_deallocate(ptr))
        {
            if (use_aligned_dealloc)
            {
                if (!use_dynamic_alignment)
                {
                    // Constant alignment
                    ::operator delete[](ptr, my_align());
                }
                else
                {
                    // Alignment depends on the number of elements
                    ::operator delete[](ptr, my_align(n));
                }
            }
            else
            {
                // Plain new
                delete[] ptr;
            }
        }
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
        "min-size",
        "INT",
        "Min size for memory pool (default: 5)"
    );
    argList::addOption
    (
        "count",
        "INT",
        "Number of elements to test (default: 10)"
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
    argList::addBoolOption
    (
        "const-align",
        "Disable dynamic alignment (default: false)"
    );

    #include "setRootCase.H"

    label count(10);

    args.readIfPresent("count", count);
    args.readIfPresent("min-size", min_size);

    use_dynamic_alignment = !args.found("const-align");
    use_aligned_alloc = !args.found("no-align-alloc");
    use_aligned_dealloc = !args.found("no-align-dealloc");


    Info<< "constant align: " << int(my_align()) << nl
        << "dynamic  align: " << int(my_align(count))
        << " for " << count << " elements" << nl
        << nl;

    {
        using T = double;
        label len = count;

        UList<T> list(my_allocate<T>(len), len);

        Info<< "List ptr: " << Foam::name(list.data()) << nl;

        list = 100;

        Info<< "List: " << list << nl;

        my_deallocate(list.data(), len);

        list = UList<T>();
    }

    return 0;
}


// ************************************************************************* //
