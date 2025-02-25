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

Application
    Test-parallel-barrier1

Description
    Simple test of local barriers communication

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "clockTime.H"
#include "IPstream.H"
#include "OPstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();
    argList::addVerboseOption();
    argList::addOption("delay", "sec", "Seconds to sleep (default 2)");

    #include "setRootCase.H"

    if (!UPstream::parRun())
    {
        Info<< "###############" << nl
            << "Not running in parallel. Stopping now" << nl
            << "###############" << endl;
        return 1;
    }

    const auto delay = args.getOrDefault<label>("delay", 2);

    Info<< nl
        << "Testing local barrier, sleep=" << delay << endl;


    const auto myProci = UPstream::myProcNo(UPstream::worldComm);
    const auto numProc = UPstream::nProcs(UPstream::worldComm);

    constexpr int uniqTag = 1516;

    clockTime timing;

    if (UPstream::master(UPstream::worldComm))
    {
        // Wait for the last rank
        UPstream::wait_done(numProc-1, UPstream::worldComm);

        // Wait for any other rank
        if (numProc > 2)
        {
            int from = UPstream::wait_done(-1, UPstream::worldComm, uniqTag);
            Pout<< "done signal from: " << from << endl;
        }
    }
    else if (myProci == numProc-1)
    {
        Foam::sleep(delay);
        UPstream::send_done(UPstream::masterNo(), UPstream::worldComm);
    }

    // Cascade sequencing (and delays)
    if (numProc > 7)
    {
        if (myProci == 2)
        {
            Foam::sleep(2*delay);
            UPstream::send_done(4, UPstream::worldComm);
        }
        else if (myProci == 4)
        {
            UPstream::wait_done(2, UPstream::worldComm);
            Foam::sleep(2*delay);
            UPstream::send_done(5, UPstream::worldComm);
        }
        else if (myProci == 5)
        {
            UPstream::wait_done(4, UPstream::worldComm);
        }
    }

    // Some arbitrary signaling rank
    if ((numProc > 2) && (myProci == numProc/2))
    {
        Pout<< "send done signal " << myProci << " -> 0" << endl;
        UPstream::send_done(UPstream::masterNo(), UPstream::worldComm, uniqTag);
    }

    Pout<< "done: " << timing.elapsedTime() << " s" << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
