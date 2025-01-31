/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "PstreamGlobals.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::DynamicList<bool> Foam::PstreamGlobals::pendingMPIFree_;
Foam::DynamicList<MPI_Comm> Foam::PstreamGlobals::MPICommunicators_;
Foam::DynamicList<MPI_Request> Foam::PstreamGlobals::outstandingRequests_;


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::PstreamGlobals::initCommunicator(const label index)
{
    if (FOAM_UNLIKELY(index < 0 || index > MPICommunicators_.size()))
    {
        FatalErrorInFunction
            << "PstreamGlobals out of sync with UPstream data. Problem."
            << Foam::abort(FatalError);
    }
    else if (index == MPICommunicators_.size())
    {
        // Extend storage with null values
        pendingMPIFree_.emplace_back(false);
        MPICommunicators_.emplace_back(MPI_COMM_NULL);
    }
    else
    {
        // Init with null values
        pendingMPIFree_[index] = false;
        MPICommunicators_[index] = MPI_COMM_NULL;
    }
}


// ************************************************************************* //
