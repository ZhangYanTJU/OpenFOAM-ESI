/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 PCOpt/NTUA
    Copyright (C) 2022      FOSS GP
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

#include "checkPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(checkPoint, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

checkPoint::checkPoint
(
    storageParameters& storageParams
)
:
    vars_(),
    storageParams_(storageParams),
    mesh_(storageParams.mesh()),
    storageMetrics_(),
    deltaT_(0.),
    active_(false),
    placeHolder_(true),
    a_(-1),
    checkPointLevel_(-1),
    checkPointTimeIndex_(-1),
    checkPointTime_(-1.)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void checkPoint::retrieve()
{
    if (!vars_)
    {
        FatalErrorInFunction
            << "Trying to retrieve solution at t = "
            << storageParams_.mesh().time().value()
            << " from checkPoint with t = "
            << checkPointTime_ << nl
            << "but no storage has been allocated." << nl
            << "placeHolder checkPoint? " << placeHolder_ << endl
            << exit(FatalError);
    }
    vars_().decompress();
}


void checkPoint::setPlaceHolder(Switch newLevel)
{
    deltaT_ = storageParams_.mesh().time().deltaTValue();
    active_ = true;
    placeHolder_ = true;
    a_ = 0;
    checkPointTimeIndex_ = storageParams_.mesh().time().timeIndex();
    checkPointTime_ = storageParams_.mesh().time().value();
    if (newLevel)
    {
        checkPointLevel_ += 1;
    }
    else
    {
        checkPointLevel_ = 0;
    }
    // Clear the  previously allocated storage
    vars_.clear();
    storageMetrics_.clear();
}


void checkPoint::store(incompressibleVars& vs, const label a)
{
    placeHolder_ = false;
    a_ = a;
    vars_.reset(new fullIncompressibleVars(vs, storageParams_, a_));
    vars_().compress();
    storageMetrics_ = vars_().storageMetrics();
}


void checkPoint::empty()
{
    active_ = false;
    placeHolder_ = true;
    checkPointTimeIndex_ = -1;
    checkPointTime_ = -1;
    checkPointLevel_ = -1;
    a_ = -1;
    vars_.clear();
    storageMetrics_.clear();
}


void checkPoint::resetToPlaceHolder()
{
    placeHolder_ = true;
    a_ = 0;
    vars_.clear();
    storageMetrics_.clear();
}


const scalarList& checkPoint::storageMetrics() const
{
    return storageMetrics_;
}


// ************************************************************************* //

} // End namespace Foam

// ************************************************************************* //
