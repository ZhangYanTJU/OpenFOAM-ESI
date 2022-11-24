/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "hostUncollatedFileOperation.H"
#include "addToRunTimeSelectionTable.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(hostUncollatedFileOperation, 0);
    addToRunTimeSelectionTable
    (
        fileOperation,
        hostUncollatedFileOperation,
        word
    );
    addToRunTimeSelectionTable
    (
        fileOperation,
        hostUncollatedFileOperation,
        comm
    );

    // Register initialisation routine. Signals need for threaded mpi and
    // handles command line arguments
    addNamedToRunTimeSelectionTable
    (
        fileOperationInitialise,
        hostUncollatedFileOperationInitialise,
        word,
        hostUncollated
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void Foam::fileOperations::hostUncollatedFileOperation::init(bool verbose)
{
    verbose = (verbose && Foam::infoDetailLevel > 0);

    if (verbose)
    {
        DetailInfo
            << "I/O    : " << this->type() << nl;

        if (ioRanks().size())
        {
            masterUncollatedFileOperation::printRanks();
        }
    }
}


Foam::fileOperations::hostUncollatedFileOperation::hostUncollatedFileOperation
(
    bool verbose
)
:
    masterUncollatedFileOperation
    (
        UPstream::allocateCommunicator
        (
            UPstream::worldComm,
            fileOperation::getGlobalSubRanks(true)  // Host
        ),
        fileOperation::getGlobalIORanks(true),  // Host
        false,  // distributedRoots
        false   // verbose
    ),
    managedComm_(comm_)
{
    init(verbose);
}


Foam::fileOperations::hostUncollatedFileOperation::hostUncollatedFileOperation
(
    const Tuple2<label, labelList>& commAndIORanks,
    const bool distributedRoots,
    bool verbose
)
:
    masterUncollatedFileOperation
    (
        commAndIORanks,
        distributedRoots,
        false   // verbose
    ),
    managedComm_(-1)  // Externally managed
{
    init(verbose);
}


Foam::fileOperations::hostUncollatedFileOperation::
hostUncollatedFileOperation
(
    const label comm,
    const labelUList& ioRanks,
    const bool distributedRoots,
    bool verbose
)
:
    masterUncollatedFileOperation(comm, ioRanks, distributedRoots, false),
    managedComm_(comm_)
{
    init(verbose);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperations::hostUncollatedFileOperation::
~hostUncollatedFileOperation()
{
    // Wait for any outstanding file operations
    flush();

    if (myComm_ != -1 && myComm_ != UPstream::worldComm)
    {
        UPstream::freeCommunicator(comm_);
    }
}


// ************************************************************************* //
