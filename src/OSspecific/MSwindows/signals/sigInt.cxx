/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2011 Symscape
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "sigInt.H"
#include "error.H"
#include "JobInfo.H"
#include "IOstreams.H"

// File-local functions
#include "signalMacros.cxx"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::sigInt::sigActive_ = false;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sigInt::sigHandler(int)
{
    resetHandler("SIGINT", SIGINT);

    JobInfo::shutdown();        // From running -> finished
    ::raise(SIGINT);            // Throw signal (to old handler)
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigInt::sigInt()
{
    set(false);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sigInt::~sigInt()
{
    unset(false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sigInt::set(bool)
{
    if (sigActive_)
    {
        return;
    }
    sigActive_ = true;

    setHandler("SIGINT", SIGINT, sigHandler);
}


void Foam::sigInt::unset(const bool verbose)
{
    if (!sigActive_)
    {
        return;
    }
    sigActive_ = false;

    resetHandler("SIGINT", SIGINT);
}


// ************************************************************************* //
