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
    Test-gather-scatter1

Description
    Simple tests for gather/scatter

\*---------------------------------------------------------------------------*/

#include "List.H"
#include "argList.H"
#include "Time.H"
#include "Pstream.H"
#include "IOstreams.H"

using namespace Foam;


//- Ostensibly the inverse of gatherList, but actually works like
//- a broadcast that skips overwriting the local rank!
template<class T>
void real_scatterList
(
    //! [in,out]
    UList<T>& values,
    [[maybe_unused]] const int tag = UPstream::msgType(),
    const int communicator = UPstream::worldComm
)
{
    if (!UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return;
    }
    else if constexpr (is_contiguous_v<T>)
    {
        // This part is a real in-place scatter:

        // In-place scatter for contiguous types - one element per rank
        // - on master:
        //   * send pointer is the full list
        //   * recv pointer is first destination
        // - on rank:
        //   * send pointer is irrelevant
        //   * recv pointer is destination in the list
        //
        // So can simply use identical pointers for send/recv

        auto* ptr = values.data() + UPstream::myProcNo(communicator);
        UPstream::mpiScatter(ptr, ptr, 1, communicator);
    }
    else
    {
        // Communication order
        const auto& commOrder = UPstream::whichCommunication(communicator);

        Pstream::scatterList_algorithm(commOrder, values, tag, communicator);
    }
}


//- gatherList_algorithm, but with specific communication style
template<class T>
void gatherList_algo
(
    const bool linear,
    //! [in,out]
    UList<T>& values,
    [[maybe_unused]] const int tag = UPstream::msgType(),
    const int communicator = UPstream::worldComm
)
{
    if (UPstream::is_parallel(communicator))
    {
        Pstream::gatherList_algorithm
        (
            UPstream::whichCommunication(communicator, linear),
            values,
            tag,
            communicator
        );
    }
}


//- scatterList_algorithm, but with specific communication style
template<class T>
void scatterList_algo
(
    const bool linear,
    //! [in,out]
    UList<T>& values,
    [[maybe_unused]] const int tag = UPstream::msgType(),
    const int communicator = UPstream::worldComm
)
{
    if (UPstream::is_parallel(communicator))
    {
        Pstream::scatterList_algorithm
        (
            UPstream::whichCommunication(communicator, linear),
            values,
            tag,
            communicator
        );
    }
}


// Perform tests
template<class ListType, class ResetCode>
void doTest(ResetCode reset)
{
    ListType values;

    reset(values);

    Pout<< nl << "before:" << flatOutput(values) << endl;
    Pstream::broadcastList(values);
    Pout<< "broadcast:" << flatOutput(values) << endl;

    reset(values);

    Pout<< nl << "before:" << flatOutput(values) << endl;
    Pstream::scatterList(values);
    Pout<< "scatter:" << flatOutput(values) << endl;

    reset(values);

    Pout<< "before:" << flatOutput(values) << endl;
    real_scatterList(values);
    Pout<< "inplace:" << flatOutput(values) << endl;

    using control = std::pair<int, int>;

    const char* algoType[2] = { "tree", "linear" };

    for
    (
        auto [gather, scatter] :
        {
            control{0, 0},
            control{1, 1},
            control{0, 1},
            control{1, 0}
        }
    )
    {
        reset(values);

        Pout<< nl << "before:" << flatOutput(values) << endl;

        gatherList_algo(gather, values);
        Pout<< "gather[" << algoType[gather] << "]:"
            << flatOutput(values) << endl;

        scatterList_algo(scatter, values);
        Pout<< "scatter[" << algoType[scatter] << "]:"
            << flatOutput(values) << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();
    argList::addVerboseOption("increase UPstream::debug level");

    #include "setRootCase.H"

    const int optVerbose = args.verbose();

    if (optVerbose)
    {
        UPstream::debug = optVerbose;
    }

    Pout<< nl << "Test contiguous" << endl;
    {
        doTest<labelList>
        (
            [](auto& values){
                if (UPstream::master())
                {
                    values = identity(UPstream::nProcs());
                }
                else
                {
                    values.resize(UPstream::nProcs());
                    values = -1;
                    values[UPstream::myProcNo()] = 10 * UPstream::myProcNo();
                }
            }
        );
    }

    Pout<< nl << "Test non-contiguous" << endl;
    {
        doTest<wordList>
        (
            [](auto& values) {
                values.resize(UPstream::nProcs());
                if (UPstream::master())
                {
                    forAll(values, i)
                    {
                        values[i] = "proc" + Foam::name(i);
                    }
                }
                else
                {
                    values = "none";
                    values[UPstream::myProcNo()] =
                        "_" + Foam::name(UPstream::myProcNo());
                }
            }
        );
    }

    // Test dummy broadcast as well
    Pout<< nl << "Test broadcastList" << endl;
    {
        wordList list;

        Pout<< nl << "before: " << flatOutput(list) << endl;

        Pstream::broadcastList(list);
        Pout<< "-> " << flatOutput(list) << endl;
    }

    // Test in-place reduce
    Pout<< nl << "Test in-place reduce" << endl;
    {
        FixedList<label, 6> list;
        list = UPstream::myProcNo();

        Pout<< nl << "before: " << flatOutput(list) << endl;

        UPstream::mpiReduce
        (
            list.data(),
            list.size(),
            UPstream::opCodes::op_sum,
            UPstream::worldComm
        );

        Pout<< "-> " << flatOutput(list) << endl;
    }

    Info<< nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
