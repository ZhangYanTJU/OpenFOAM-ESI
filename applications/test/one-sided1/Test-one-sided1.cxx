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
    Test-one-sided1

Description
    Simple test of one-sided communication

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IPstream.H"
#include "OPstream.H"
#include "SubField.H"
#include "vector.H"
#include "IOstreams.H"

using namespace Foam;

template<class T>
Ostream& printSpanInfo(Ostream& os, const UList<T>& span)
{
    os  << "addr=" << Foam::name(span.cdata())
        << " size= " << span.size();

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();
    argList::addVerboseOption();
    argList::addBoolOption("no-shared", "disable shared memory tests");
    argList::addBoolOption("no-sleep", "disable sleep for async test");

    #include "setRootCase.H"

    const bool with_shared = !args.found("no-shared");
    const bool with_sleep = !args.found("no-sleep");

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl
        << "nProcs = " << UPstream::nProcs()
        << " with " << UPstream::nComms() << " predefined comm(s)" << nl;

    if (!UPstream::parRun())
    {
        Info<< "###############" << nl
            << "Not running in parallel. Stopping now" << nl
            << "###############" << endl;
        return 1;
    }

    const auto myProci = UPstream::myProcNo();
    const auto numProc = UPstream::nProcs();

    // Make some windows
    Field<label> buffer(10 + myProci);
    buffer = myProci;

    Pout<< "input: " << flatOutput(buffer) << endl;

    UPstream::Window win;
    win.create(buffer, UPstream::worldComm);

    // Pass 1
    // - grab things from sub-ranks
    if (UPstream::master())
    {
        win.lock_all(true);

        win.get
        (
            buffer.slice(4, 2),
            1,  // target_rank
            2   // target_disp
        );

        win.unlock_all();
    }

    Pout<< "output: " << flatOutput(buffer) << endl;

    // Pass 2:
    // accumulate into master
    if (UPstream::is_subrank())
    {
        win.lock(0);

        win.put
        (
            UPstream::opCodes::op_sum,
            buffer.slice(2, 4),
            UPstream::masterNo(),
            2  // target_disp
        );

        win.unlock(0);
    }

    Pout<< "updated: " << flatOutput(buffer) << endl;

    // Pass 3:
    // Update some values - something very asynchronous
    if (UPstream::is_subrank())
    {
        if (with_sleep)
        {
            if (UPstream::myProcNo() % 3)
            {
                Foam::sleep(3);
            }
            else
            {
                Foam::sleep(1);
            }
        }
        buffer *= 10;
        forAll(buffer, i)
        {
            buffer[i] *= 1 + (i % 3);
        }
    }

    // Needs a process sync, otherwise master fetches old values
    UPstream::barrier(UPstream::worldComm);

    label lastValue(-1);

    if (UPstream::master())
    {
        win.lock_all(true);

        for (const auto proci : UPstream::subProcs())
        {
            win.fetch_and_op
            (
                UPstream::opCodes::op_sum,
                buffer[0],
                lastValue,
                proci,
                2  // target_disp
            );
        }

        // Force changes to occur
        win.flush_all();
        win.unlock_all();
    }


    Pout<< "last-value : " << lastValue << nl
        << "final : " << flatOutput(buffer) << endl;

    labelList allUpdates;

    if (UPstream::master())
    {
        allUpdates.resize(UPstream::nProcs(), -10);

        win.lock_all(true);

        for (const auto proci : UPstream::subProcs())
        {
            win.get_value
            (
                allUpdates[proci],
                proci,
                2  // target_disp
            );
        }

        win.flush_all();
        win.unlock_all();
    }

    Info<< "gets: " << flatOutput(allUpdates) << endl;


    // This should fail (runtime)
    #if 0
    if (UPstream::master())
    {
        labelPair value1(-1, -1);

        win.lock_all(true);

        for (const auto proci : UPstream::subProcs())
        {
            win.fetch_and_op
            (
                UPstream::opCodes::op_sum,
                value1,
                lastValue,
                proci,
                8  // target_disp
            );
        }

        win.unlock_all();
    }
    #endif

    // Last thing before closing out
    // replace values. Not very efficient...

    // Persistent data to move onto target:
    const label newValue(333);
    const label multiplier(-3);

    if (UPstream::master())
    {
        win.lock_all(true);

        for (const auto proci : UPstream::subProcs())
        {
            win.fetch_and_op
            (
                UPstream::opCodes::op_replace,
                newValue,
                lastValue,
                proci,  // target_rank
                3       // target_disp
            );

            win.put_value
            (
                UPstream::opCodes::op_prod,
                multiplier,
                proci,  // target_rank
                5       // target_disp
            );
        }
        win.unlock_all();
    }

    win.close();  // process collective

    Pout<< "modified: " << flatOutput(buffer) << endl;

    if (with_shared)
    {
        // Make some shared window
        UList<label> newBuffer;

        {
            label localLen(0);

            if
            (
                (myProci == 3)
             || (myProci == numProc-2)
            )
            {
                localLen = 0;
            }
            else
            {
                localLen = (10 + UPstream::myProcNo());
            }

            // Just to prove that we can shallow copy the view...
            newBuffer =
                win.allocate_shared<label>(localLen, UPstream::worldComm);
        }

        newBuffer = UPstream::myProcNo();

        Pout<< "Shared: " << flatOutput(newBuffer) << endl;

        {
            UList<label> local = win.view<label>();
            Pout<< "local win: ";
            printSpanInfo(Pout, local) << endl;
        }

        Pout<< "Query rank1" << endl;

        {
            // UPtrList<UList<label>> totalList(UPstream::nProcs());
            //
            // totalList.set(0, &newBuffer);

            const label* ptr0 = nullptr;

            {
                UList<label> buf = win.view_shared<label>(0);
                ptr0 = buf.cdata();

                Pout<< "addr 0 = " << Foam::name(ptr0)
                    << " diff = " << label(0)
                    << " + " << buf.size() << endl;
            }

            // UList<label> other = win.global_view<label>();

            for (const auto proci : UPstream::subProcs())
            {
                UList<label> other = win.view_shared<label>(proci);

                const label* ptr = other.cdata();

                Pout<< "addr " << proci << " = "
                    << Foam::name(ptr)
                    << " diff = " << label(ptr - ptr0)
                    << " + " << other.size() << endl;
                // totalList.set(proci, &other);
            }
        }

        win.close();
    }

    // Since close() is ignored on null window,
    // can call it an arbitrary number of times
    win.close();
    win.close();
    win.close();

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
