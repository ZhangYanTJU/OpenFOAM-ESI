/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2016-2025 OpenCFD Ltd.
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

#include "memInfo.H"
#include "IOstreams.H"
// #include "OSspecific.H"  // For pid()

#undef DebugInfo        // Windows name clash with OpenFOAM messageStream

#include <cstdlib>
#include <fstream>
#include <string>

#define WIN32_LEAN_AND_MEAN
#ifdef FOAM_USE_WINDOWS_PSAPI
#include <windows.h>
#include <psapi.h>
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::memInfo::memInfo()
:
    peak_(0),
    size_(0),
    rss_(0),
    free_(0)
{
    populate();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::memInfo::good() const noexcept
{
    return peak_ > 0;
}


void Foam::memInfo::clear() noexcept
{
    peak_ = size_ = rss_ = free_ = 0;
}


void Foam::memInfo::populate()
{
    #ifdef FOAM_USE_WINDOWS_PSAPI
    HANDLE proc = ::OpenProcess
    (
        PROCESS_QUERY_INFORMATION | PROCESS_VM_READ,
        FALSE,
        ::GetCurrentProcessId()  // This is Foam::pid()
    );

    if (proc)
    {
        if
        (
            PROCESS_MEMORY_COUNTERS pmc;
            ::GetProcessMemoryInfo(proc, &pmc, sizeof(pmc))
        )
        {
            // Convert from bytes -> kibi-bytes (1024)
            peak_ = int64_t(pmc.PeakWorkingSetSize) / 1024;
            size_ = int64_t(pmc.WorkingSetSize) / 1024;
            rss_ = int64_t(pmc.WorkingSetSize - pmc.PagefileUsage) / 1024;
        }
        CloseHandle(proc);

        if
        (
            PERFORMANCE_INFORMATION pinfo;
            ::GetPerformanceInfo(&pinfo, sizeof(pinfo))
        )
        {
            // Convert from pages -> bytes -> kibi-bytes (1024)
            free_ = int64_t(pinfo.PhysicalAvailable*pinfo.PageSize) / 1024;
        }
    }
    #endif  /* FOAM_USE_WINDOWS_PSAPI */
}


const Foam::memInfo& Foam::memInfo::update()
{
    clear();
    populate();
    return *this;
}


void Foam::memInfo::writeEntries(Ostream& os) const
{
    os.writeEntry("size", size_);
    os.writeEntry("peak", peak_);
    os.writeEntry("rss", rss_);
    os.writeEntry("free", free_);
    os.writeEntry("units", "kB");  // kibi-bytes (1024)
}


void Foam::memInfo::writeEntry(const word& keyword, Ostream& os) const
{
    os.beginBlock(keyword);
    writeEntries(os);
    os.endBlock();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// Foam::Istream& Foam::operator>>(Istream& is, memInfo& m)
// {
//     is.readBegin("memInfo");
//     is  >> m.peak_ >> m.size_ >> m.rss_ >> m.free_;
//     is.readEnd("memInfo");
//
//     is.check(FUNCTION_NAME);
//     return is;
// }


Foam::Ostream& Foam::operator<<(Ostream& os, const memInfo& m)
{
    os  << token::BEGIN_LIST
        << m.peak() << token::SPACE
        << m.size() << token::SPACE
        << m.rss()  << token::SPACE
        << m.free()
        << token::END_LIST;

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
