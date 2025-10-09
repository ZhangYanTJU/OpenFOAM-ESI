/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2025 OpenCFD Ltd.
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

Application
    foamListRegions

Group
    grpPostProcessingUtilities

Description
    List volume regions from constant/regionProperties
    or area regions from constant/finite-area/regionProperties

Usage
    \b foamListRegions [OPTION]

Note
    The OpenFOAM banner information is suppressed so that the output can be
    piped into another command.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "regionProperties.H"

// Same as faMesh::prefix() but without additional linkage
constexpr const char* const faMeshPrefix = "finite-area";

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "List volume regions from constant/regionProperties,\n"
        "or area regions from constant/finite-area/regionProperties"
    );

    argList::noBanner();
    argList::noParallel();
    argList::noJobInfo();
    argList::noFunctionObjects();  // Never use function objects
    // No profiling since there is no time loop

    argList::addBoolOption
    (
        "finite-area",
        "List constant/finite-area/regionProperties (if available)"
    );

    argList::addBoolOption
    (
        "optional",
        "A missing regionProperties is not treated as an error"
    );

    argList::addDryRunOption
    (
        "Make reading optional and add verbosity"
    );
    argList::addVerboseOption("Additional verbosity");

    // Arguments are optional (non-mandatory)
    argList::noMandatoryArgs();
    argList::addArgument("regionType ... regionType");

    #include "setRootCase.H"

    const bool dryRun = args.dryRun();
    int optVerbose = args.verbose();

    if (dryRun && !optVerbose)
    {
        ++optVerbose;
    }

    // File is optional, not an error
    const bool isOptional = args.found("optional");

    // Use finite-area regions
    const bool doFiniteArea = args.found("finite-area");

    // The number of optional region filters to apply
    const label nFilters = (args.size()-1);

    IOobjectOption::readOption readOpt(IOobjectOption::MUST_READ);

    if (dryRun || isOptional || doFiniteArea)
    {
        // The finite-area regionProperties are also considered optional
        readOpt = IOobjectOption::READ_IF_PRESENT;
    }

    // Silent version of "createTime.H", without libraries
    Time runTime
    (
        Time::controlDictName,
        args,
        false,  // no enableFunctionObjects
        false   // no enableLibs
    );

    regionProperties regionProps;
    if (doFiniteArea)
    {
        regionProps = regionProperties(runTime, faMeshPrefix, readOpt);
    }
    else
    {
        regionProps = regionProperties(runTime, readOpt);
    }

    // Some reporting...
    if (regionProps.empty())
    {
        if (doFiniteArea)
        {
            InfoErr<< "No finite-area region types" << nl;
        }
        else if (isOptional)
        {
            InfoErr<< "No region types" << nl;
        }
    }
    else if (optVerbose)
    {
        InfoErr << "Have " << regionProps.size();

        if (doFiniteArea)
        {
            InfoErr<< " finite-area";
        }
        InfoErr
            << " region types, "
            << regionProps.count() << " regions" << nl << nl;
    }


    // We now handle checking args and general sanity etc.

    DynamicList<word> regionTypes;

    if (isOptional && regionProps.empty())
    {
        // Nothing to do...
    }
    else if (nFilters > 0)
    {
        // Apply region filters

        regionTypes.reserve_exact
        (
            Foam::min(nFilters, regionProps.size())
        );

        // No duplicates, and no duplicate warnings
        wordHashSet uniq;

        for (label argi = 1; argi < args.size(); ++argi)
        {
            word regType(args[argi]);

            if (uniq.insert(regType))
            {
                if (regionProps.contains(regType))
                {
                    if (!regionTypes.contains(regType))
                    {
                        regionTypes.push_back(std::move(regType));
                    }
                }
                else
                {
                    InfoErr<< "No region-type: " << regType << nl;
                }
            }
        }
    }
    else
    {
        // Take all regions
        regionTypes = regionProps.sortedToc();
    }


    for (const word& regionType : regionTypes)
    {
        if (const auto iter = regionProps.cfind(regionType); iter.good())
        {
            for (const word& regionName : iter.val())
            {
                Info<< regionName << nl;
            }
        }
    }

    return 0;
}


// ************************************************************************* //
