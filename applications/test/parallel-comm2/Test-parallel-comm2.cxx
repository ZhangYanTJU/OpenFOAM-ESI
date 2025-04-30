/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2025 OpenCFD Ltd.
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
    Test-parallel-comm2

Description
    Basic communicator tests

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IPstream.H"
#include "OPstream.H"
#include "Pair.H"
#include "Tuple2.H"
#include "IOstreams.H"
#include "PstreamReduceOps.H"
#include "SHA1.H"

#include "openfoam_mpi.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void rankInfo(const label comm)
{
    const int ranki = UPstream::myProcNo(comm);

    Pout<< "comm:" << comm
        << "(parent:" << UPstream::parent(comm) << ')'
        << " rank:" << ranki
        << "(sub:" << UPstream::is_subrank(comm)
        << ") nProcs:" << UPstream::nProcs(comm);
      // << " baseProcNo:" << UPstream::baseProcNo(comm, ranki);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noCheckProcessorDirectories();
    argList::addBoolOption("info", "information");
    argList::addBoolOption("print-tree", "Report tree(s) as graph");
    argList::addBoolOption("no-test", "Disable general tests");
    argList::addBoolOption("split", "Test Pstream split-comm");
    argList::addBoolOption("host-comm", "Test Pstream host-comm");
    argList::addBoolOption("host-broadcast", "Test host-base broadcasts");

    #include "setRootCase.H"

    const bool optPrintTree = args.found("print-tree");
    bool generalTest = !args.found("no-test");

    Info<< nl
        << "parallel:" << UPstream::parRun()
        << " nProcs = " << UPstream::nProcs()
        << " with " << UPstream::nComms() << " predefined comm(s)."
        << " proc:" << UPstream::myProcNo() << nl;

    if (UPstream::parRun() && optPrintTree)
    {
        // Info<< "comms: "
        //     << UPstream::whichCommunication(UPstream::worldComm) << nl;
        UPstream::printCommTree(UPstream::commWorld());
    }

    if (UPstream::parRun())
    {
        Pout<< "world ranks: 0.."
            << UPstream::nProcs(UPstream::commWorld())-1 << nl;

        Pout<< "inter-node ranks: " << UPstream::numNodes() << ' '
            << flatOutput(UPstream::procID(UPstream::commInterNode())) << nl;

        Pout<< "local-node ranks: "
            << flatOutput(UPstream::procID(UPstream::commLocalNode())) << nl;
    }

    if (UPstream::parRun() && args.found("split"))
    {
        Info<< "split: alternative ranks" << nl;

        const auto myRank = UPstream::myProcNo();

        int colour =
        (
            (myRank == 5 || myRank == 6)  // Exclude these ones
          ? -1
          : (myRank % 2)
        );

        UPstream::communicator comm =
            UPstream::communicator::split(UPstream::commWorld(), colour, true);

        Pout<< "split ranks (colour=" << colour << ") "
            << flatOutput(UPstream::procID(comm.comm())) << nl;

        comm.reset();
        comm =
            UPstream::communicator::split(UPstream::commWorld(), colour, false);

        Pout<< "Split ranks (colour=" << colour << ") "
            << flatOutput(UPstream::procID(comm.comm())) << nl;
    }


    if (args.found("info"))
    {
        Info<< nl;

        // Process IDs within a given communicator
        for (label comm = 0; comm < UPstream::nComms(); ++comm)
        {
            Info<< "comm=" << comm
                << " procIDs: " << flatOutput(UPstream::procID(comm)) << endl;
            rankInfo(comm);
            Pout<< nl;
        }
        Pout<< endl;
    }

    if (UPstream::parRun() && args.found("host-comm"))
    {
        generalTest = false;
        Info<< nl << "[pstream host-comm]" << nl << endl;

        const label commInterNode = UPstream::commInterNode();
        const label commLocalNode = UPstream::commLocalNode();

        Pout<< "Host rank " << UPstream::myProcNo(commLocalNode)
            << " / " << UPstream::nProcs(commLocalNode)
            << " on " << hostName()
            << ", inter-rank: " << UPstream::myProcNo(commInterNode)
            << " / " << UPstream::nProcs(commInterNode)
            << ", host leader:" << UPstream::master(commInterNode)
            << " sub-rank:" << UPstream::is_subrank(commInterNode)
            << endl;

        {
            // Info<< "host-master: "
            //     << UPstream::whichCommunication(commInterNode) << endl;

            UPstream::printCommTree(commInterNode);
            UPstream::printCommTree(commLocalNode);
        }
    }

    if (UPstream::parRun() && args.found("host-broadcast"))
    {
        generalTest = false;
        Info<< nl << "[pstream host-broadcast]" << nl << endl;

        const label commInterNode = UPstream::commInterNode();
        const label commLocalNode = UPstream::commLocalNode();

        Pout<< "world rank: " << UPstream::myProcNo(UPstream::commWorld())
            << " host-leader rank: "
            << UPstream::myProcNo(UPstream::commInterNode())
            << " intra-host rank: "
            << UPstream::myProcNo(UPstream::commLocalNode())
            << endl;

        label value1(0), value2(0), value3(0);
        label hostIndex = UPstream::myProcNo(commInterNode);

        if (UPstream::master(commInterNode))
        {
            value1 = 100;
            value2 = 200;
        }
        if (UPstream::master(commLocalNode))
        {
            value3 = 300;
        }

        Pstream::broadcast(value1, commInterNode);
        Pstream::broadcast(value2, commLocalNode);
        Pstream::broadcast(hostIndex, commLocalNode);

        Pout<< "host: " << hostIndex
            << " broadcast 1: "
            << value1 << ' '
            << value2 << ' '
            << value3 << endl;

        // re-broadcast
        Pstream::broadcast(value1, commLocalNode);
        Pout<< "host: " << hostIndex
            << " broadcast 2: "
            << value1 << endl;


        label reduced1 = value1;
        label reduced2 = value1;

        Foam::reduce
        (
            reduced1,
            sumOp<label>(),
            UPstream::msgType(),
            commLocalNode
        );

        Foam::reduce
        (
            reduced2,
            sumOp<label>(),
            UPstream::msgType(),
            commInterNode
        );

        Pout<< "value1: (host) " << reduced1
            << " (leader) " << reduced2 << endl;

        // Pout<< "ranks: " << UPstream::nProcs(commInterNode) << endl;

        wordList strings;
        if (UPstream::is_rank(commInterNode))
        {
            strings.resize(UPstream::nProcs(commInterNode));
            strings[UPstream::myProcNo(commInterNode)] = name(pid());
        }

        // Some basic gather/scatter
        Pstream::allGatherList(strings, UPstream::msgType(), commInterNode);

        Pout<< "pids " << flatOutput(strings) << endl;

        Foam::reverse(strings);

        Pstream::broadcast(strings, commLocalNode);
        Pout<< "PIDS " << flatOutput(strings) << endl;
    }


    if (UPstream::parRun() && generalTest)
    {
        #if 1
        // With first ranks
        labelList subRanks =
            identity(UPstream::nProcs(UPstream::commWorld()) / 2);

        UPstream::communicator newComm;

        newComm.reset(UPstream::commWorld(), subRanks);
        label localRanki = UPstream::myProcNo(newComm);

        const int myProci = UPstream::myProcNo(UPstream::commWorld());

        Pout.prefix() =
        (
            '[' + Foam::name(myProci) + " a:" + Foam::name(localRanki) + "] "
        );

        Info<< "procIDs: "
            << flatOutput(UPstream::procID(newComm)) << endl;

        rankInfo(newComm);
        Pout<< endl;
        #endif

        #if 1
        // With every other rank
        subRanks = identity(UPstream::nProcs(UPstream::commWorld()));

        for (label& val : subRanks)
        {
            if (val % 2) val = -1;
        }

        newComm.reset(UPstream::commWorld(), subRanks);
        localRanki = UPstream::myProcNo(newComm);

        Pout.prefix() =
        (
            '[' + Foam::name(myProci) + " b:" + Foam::name(localRanki) + "] "
        );

        Info<< "procIDs: "
            << flatOutput(UPstream::procID(newComm)) << endl;

        rankInfo(newComm);
        Pout<< endl;
        #endif
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
