/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "function1Base.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::function1Base::function1Base(const word& entryName)
:
    refCount(),
    name_(entryName),
    minLimit_(-GREAT),
    maxLimit_(GREAT)
{}


Foam::function1Base::function1Base
(
    const word& entryName,
    const dictionary& dict
)
:
    refCount(),
    name_(entryName),
    minLimit_(dict.getOrDefault<scalar>("minLimit", -GREAT)),
    maxLimit_(dict.getOrDefault<scalar>("maxLimit", GREAT))
{}


Foam::function1Base::function1Base(const function1Base& rhs)
:
    refCount(),
    name_(rhs.name_),
    minLimit_(rhs.minLimit_),
    maxLimit_(rhs.maxLimit_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::function1Base::convertTimeBase(const Time& t)
{
    minLimit_ = t.userTimeToTime(minLimit_);
    maxLimit_ = t.userTimeToTime(maxLimit_);
}


Foam::scalar Foam::function1Base::limitValue(const scalar x) const
{
    return min(max(minLimit_, x), maxLimit_);
}


// TBD
// void Foam::function1Base::writeEntries(Ostream& os) const
// {
//     if (minLimit_ > -ROOTGREAT)
//     {
//         os.writeEntry("minLimit", minLimit_);
//     }
//     if (maxLimit_ < ROOTGREAT)
//     {
//         os.writeEntry("maxLimit", maxLimit_);
//     }
// }


// ************************************************************************* //
