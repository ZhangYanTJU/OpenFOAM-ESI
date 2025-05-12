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

#include "MemoryPool.H"
#include "debug.H"
#include "dictionary.H"
#include "sigFpe.H"
#include "OSspecific.H"  // For getEnv

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#ifdef FOAM_USE_UMPIRE

// #include <cerrno>
#include <cinttypes>
#include <tuple>

#include "umpire/Allocator.hpp"
#include "umpire/ResourceManager.hpp"
#include "umpire/strategy/AlignedAllocator.hpp"
#include "umpire/strategy/DynamicPoolList.hpp"

static bool disabled_(false);
static umpire::Allocator aligned_allocator;
static umpire::Allocator pooled_allocator;
static umpire::ResourceManager* manager_(nullptr);
static umpire::ResourceManager* suspended_(nullptr);

#endif

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

#ifdef FOAM_USE_UMPIRE

namespace
{

// Different supported allocation types
enum class Types { undefined, none, host, device, managed };

typedef std::tuple<Types, std::size_t, std::size_t> ctrlTuple;

// Extract key=INT, the key includes the '='
int getIntParameter(const std::string& key, const std::string& ctrl)
{
    int val(0);

    const auto pos = ctrl.find(key);

    if (pos == std::string::npos)
    {
        return val;
    }

    const char* buf = (ctrl.data() + pos + key.size());

    char *endptr = nullptr;
    errno = 0;
    auto parsed = std::strtoimax(buf, &endptr, 10);

    if (errno || endptr == buf)
    {
        // Some type of error OR no conversion
    }
    else
    {
        val = int(parsed);
    }

    return val;
}


ctrlTuple getControlValues(const std::string& ctrl)
{
    ctrlTuple result(Types::undefined, 0, 0);

    bool checkParam = false;

    // Also find things that look like Switch constants.
    // Unfortunately need to do this manually since Switch::find()
    // itself would not manage to parse something like "true; size=10"

    if (ctrl.empty())
    {
        // Nothing => undefined
    }
    else if
    (
        std::string::npos != ctrl.find("false")     // ctrl.contains("false")
     || std::string::npos != ctrl.find("off")       // ctrl.contains("off")
     || std::string::npos != ctrl.find("no")        // ctrl.contains("no")
     || std::string::npos != ctrl.find("none")      // ctrl.contains("none")
    )
    {
        std::get<0>(result) = Types::none;
    }
    else if
    (
        std::string::npos != ctrl.find("true")      // ctrl.contains("true")
     || std::string::npos != ctrl.find("on")        // ctrl.contains("on")
     || std::string::npos != ctrl.find("yes")       // ctrl.contains("yes")

     || std::string::npos != ctrl.find("host")      // ctrl.contains("host")
     || std::string::npos != ctrl.find("system")    // ctrl.contains("system")
    )
    {
        std::get<0>(result) = Types::host;
        checkParam = true;
    }

    // These need more testing
    else if
    (
        std::string::npos != ctrl.find("device")    // ctrl.contains("device")
    )
    {
        std::get<0>(result) = Types::device;
        checkParam = true;
    }
    else if
    (
        std::string::npos != ctrl.find("managed")   // ctrl.contains("managed")
    )
    {
        std::get<0>(result) = Types::managed;
        checkParam = true;
    }

    if (checkParam)
    {
        std::get<1>(result) = getIntParameter("size=", ctrl);
        std::get<2>(result) = getIntParameter("incr=", ctrl);
    }

    return result;
}


bool create_from(const ctrlTuple& controls, bool verbose)
{
    using namespace Foam;

    if (manager_ || suspended_)
    {
        // Already created
        return true;
    }

    // Type, initial size, increment
    auto [which, size, incr] = controls;

    // std::cerr
    //     << "which=" << int(which)
    //     << ", size=" << int(size)
    //     << ", incr=" << int(incr) << '\n';


    constexpr size_t MegaByte(1024*1024);

    switch (which)
    {
        case Types::undefined :
        {
            if (verbose)
            {
                Info<< "memory pool : unused" << nl;
            }
            break;
        }

        case Types::none :
        {
            if (verbose)
            {
                Info<< "memory pool : disabled" << nl;
            }
            break;
        }

        case Types::host :
        {
            // Default sizing parameters
            if (!size) size = 1024;
            if (!incr) incr = 5;

            auto& rm = umpire::ResourceManager::getInstance();
            manager_ = &rm;

            aligned_allocator =
                rm.makeAllocator<umpire::strategy::AlignedAllocator>
                (
                    "aligned_allocator",
                    rm.getAllocator("HOST"),

                    // alignment
                    256
                );

            pooled_allocator =
                rm.makeAllocator<umpire::strategy::DynamicPoolList>
                (
                    "openfoam_HOST_pool",
                    aligned_allocator,

                    // initial block allocation size
                    (size*MegaByte),

                    // incremental block allocation size
                    (incr*MegaByte)
                );

            if (verbose)
            {
                Info<< "memory pool : host (size="
                    << int(size) << "MB, incr="
                    << int(incr) << "MB)\n";
            }
            break;
        }

        case Types::device :
        {
            auto& rm = umpire::ResourceManager::getInstance();
            manager_ = &rm;

            aligned_allocator = rm.getAllocator("DEVICE");

            pooled_allocator =
                rm.makeAllocator<umpire::strategy::DynamicPoolList>
                (
                    "openfoam_DEVICE_pool",
                    aligned_allocator
                );

            if (verbose)
            {
                Info<< "memory pool : device" << nl;
            }
            break;
        }

        case Types::managed :
        {
            // Default sizing parameters
            if (!size) size = 10*1024;
            if (!incr) incr = 10;

            auto& rm = umpire::ResourceManager::getInstance();
            manager_ = &rm;

            aligned_allocator = rm.getAllocator("UM");

            pooled_allocator =
                rm.makeAllocator<umpire::strategy::DynamicPoolList>
                (
                    "openfoam_UM_pool",
                    aligned_allocator,

                    // initial block allocation size
                    (size*MegaByte),

                    // incremental block allocation size
                    (incr*MegaByte)
                );

            if (verbose)
            {
                Info<< "memory pool : managed (size="
                    << int(size) << "MB, incr="
                    << int(incr) << "MB)\n";
            }
            break;
        }
    }

    return (which != Types::undefined && which != Types::none);
}


} // End anonymous namespace

#endif  // FOAM_USE_UMPIRE


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

// bool Foam::MemoryPool::create(const std::string& ctrl, bool verbose)
// {
//     #ifdef FOAM_USE_UMPIRE
//     if (manager_ || suspended_)
//     {
//         // Already created
//         return true;
//     }
//
//     auto controls = getControlValues(ctrl);
//
//     return create_from(controls, verbose);
//     #else
//     return false;
//     #endif
// }


bool Foam::MemoryPool::create(bool verbose)
{
    #ifdef FOAM_USE_UMPIRE
    if (disabled_)
    {
        // Disallowed
        return false;
    }
    else if (manager_ || suspended_)
    {
        // Already created
        return true;
    }

    // First check environment
    auto controls = getControlValues(Foam::getEnv("FOAM_MEMORY_POOL"));

    if (std::get<0>(controls) == Types::none)
    {
        // Disabled from environment - has highest priority
        disabled_ = true;
    }

    // Currently no easy way to handle <system>/controlDict...

    // Fallback from etc/controlDict
    if (std::get<0>(controls) == Types::undefined)
    {
        // From central etc/controlDict
        const auto& dict = Foam::debug::optimisationSwitches();

        if (auto* eptr = dict.findStream("memory_pool", keyType::LITERAL))
        {
            const token& firstToken = eptr->front();

            if (firstToken.isStringType())
            {
                controls = getControlValues(firstToken.stringToken());
            }
        }
    }

    return create_from(controls, verbose);
    #else
    if (verbose)
    {
        Info<< "memory pool : not available" << nl;
    }
    return false;
    #endif
}


void Foam::MemoryPool::destroy(bool verbose)
{
    // Nothing currently needed but could add in something like this:

    // if (manager_ || suspended_)
    // {
    //     pooled_allocator.release();
    // }

    // However, need to find the proper sequence within
    // Foam::exit() or UPstream::exit() ...
}


bool Foam::MemoryPool::active() noexcept
{
    #ifdef FOAM_USE_UMPIRE
    return bool(manager_);
    #else
    return false;
    #endif
}


bool Foam::MemoryPool::suspend() noexcept
{
    #ifdef FOAM_USE_UMPIRE
    bool status(suspended_);
    if (manager_)  // <- and (!suspended_)
    {
        std::swap(manager_, suspended_);
    }
    return status;
    #else
    return false;
    #endif
}


void Foam::MemoryPool::resume() noexcept
{
    #ifdef FOAM_USE_UMPIRE
    if (suspended_)  // <- and (!manager_)
    {
        std::swap(manager_, suspended_);
    }
    #endif
}


bool Foam::MemoryPool::is_pool(void* ptr)
{
    #ifdef FOAM_USE_UMPIRE
    if (ptr)
    {
        if (manager_)
        {
            return manager_->hasAllocator(ptr);
        }
        else if (suspended_)
        {
            return suspended_->hasAllocator(ptr);
        }
    }
    #endif

    return false;
}


void* Foam::MemoryPool::try_allocate(std::size_t nbytes)
{
    void* ptr = nullptr;

    #ifdef FOAM_USE_UMPIRE
    if (manager_)
    {
        ptr = pooled_allocator.allocate(nbytes);

        // std::cerr<< "allocate(" << int(nbytes) << ")\n";

        // Optionally fill with NaN (depends on current flags)
        Foam::sigFpe::fillNan_if(ptr, nbytes);

        if (!ptr)
        {
            // Pout<< "umpire failed to allocate memory\n";
        }
    }
    #endif

    return ptr;
}


bool Foam::MemoryPool::try_deallocate(void* ptr)
{
    #ifdef FOAM_USE_UMPIRE
    if (ptr)
    {
        if (manager_)
        {
            if (manager_->hasAllocator(ptr))  // <- ie, is_pool()
            {
                // std::cerr<< "deallocate()\n";
                manager_->deallocate(ptr);
                return true;
            }
        }
        else if (suspended_)
        {
            // Deallocate even if nominally suspended

            if (suspended_->hasAllocator(ptr))  // <- ie, is_pool()
            {
                // std::cerr<< "deallocate()\n";
                suspended_->deallocate(ptr);
                return true;
            }
        }
    }
    #endif

    return (!ptr);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
