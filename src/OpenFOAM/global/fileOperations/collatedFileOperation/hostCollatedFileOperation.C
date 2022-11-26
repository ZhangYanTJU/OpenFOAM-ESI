/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenFOAM Foundation
    Copyright (C) 2021-2022 OpenCFD Ltd.
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

#include "hostCollatedFileOperation.H"
#include "addToRunTimeSelectionTable.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(hostCollatedFileOperation, 0);
    addToRunTimeSelectionTable
    (
        fileOperation,
        hostCollatedFileOperation,
        word
    );

    // Register initialisation routine. Signals need for threaded mpi and
    // handles command line arguments
    addNamedToRunTimeSelectionTable
    (
        fileOperationInitialise,
        hostCollatedFileOperationInitialise,
        word,
        hostCollated
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void Foam::fileOperations::hostCollatedFileOperation::init(bool verbose)
{
    verbose = (verbose && Foam::infoDetailLevel > 0);

    if (verbose)
    {
        this->printBanner(ioRanks_.size());
    }
}


Foam::fileOperations::hostCollatedFileOperation::hostCollatedFileOperation
(
    bool verbose
)
:
    collatedFileOperation
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


Foam::fileOperations::hostCollatedFileOperation::hostCollatedFileOperation
(
    const Tuple2<label, labelList>& commAndIORanks,
    const bool distributedRoots,
    bool verbose
)
:
    collatedFileOperation
    (
        commAndIORanks,
        distributedRoots,
        false   // verbose
    ),
    managedComm_(-1)  // Externally managed
{
    if (verbose && Foam::infoDetailLevel > 0)
    {
        this->printBanner(ioRanks_.size());
    }
}


Foam::fileOperations::hostCollatedFileOperation::hostCollatedFileOperation
(
    const label comm,
    const labelUList& ioRanks,
    const bool distributedRoots,
    bool verbose
)
:
    collatedFileOperation
    (
        comm,
        ioRanks,
        distributedRoots,
        false  // verbose
    ),
    managedComm_(-1)  // Externally managed
{
    if (verbose && Foam::infoDetailLevel > 0)
    {
        this->printBanner(ioRanks_.size());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperations::hostCollatedFileOperation::~hostCollatedFileOperation()
{
    if (managedComm_ >= 0 && managedComm_ != UPstream::worldComm)
    {
        UPstream::freeCommunicator(managedComm_);
    }
}


// ************************************************************************* //
