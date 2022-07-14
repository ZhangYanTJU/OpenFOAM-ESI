/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 PCOpt/NTUA
    Copyright (C) 2022 FOSS GP
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

#include "error.H"
#include "unsteadyTimeManipulation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(unsteadyTimeManipulation, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::unsteadyTimeManipulation::printOut(const char* functionName)
{
    if (debug)
    {
        Info<< "After call to " << functionName << endl;
        OSstream& os = Info.stream();
        os.beginBlock();
        os<< *this;
        os.endBlock();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::unsteadyTimeManipulation::unsteadyTimeManipulation
(
    const fvMesh& mesh
)
:
    MeshObject<fvMesh, UpdateableMeshObject, unsteadyTimeManipulation>(mesh),
    time_(const_cast<Time&>(mesh.time())),
    startTime_(time_.startTime()),
    startTimeIndex_(time_.startTimeIndex()),
    endTimeIndex_(-1),
    span_(time_.endTime().value() - time_.startTime().value()),
    storedTime_(nullptr),
    storedTimeIndex_(nullptr),
    storedEndTime_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::unsteadyTimeManipulation::moveToPrimalStartTime(bool setEndTime)
{
    time_.setTime(startTime_, startTimeIndex_);
    if (setEndTime)
    {
        time_.setEndTime(startTime_.value() + span_);
    }
    printOut(FUNCTION_NAME);
}


void Foam::unsteadyTimeManipulation::moveToAdjointStartTime
(
    bool setEndTime,
    bool setTimeIndex
)
{
    label index(setTimeIndex ? endTimeIndex_ : time_.timeIndex());
    time_.setTime(startTime_.value() + span_, index);
    if (setEndTime)
    {
        // Set endTime of the adjoint solver as the start of the primal one,
        // plus a deltaT
        time_.setEndTime(startTime_.value() + time_.deltaT().value());
    }
    printOut(FUNCTION_NAME);
}


void Foam::unsteadyTimeManipulation::newOptimisationCycle()
{
    // Increment endTime by a time span
    time_.setEndTime(time_.value() + span_);

    // Adjust startTime of each optimisation loop
    startTime_ = time_;
    startTimeIndex_ = time_.timeIndex();
    printOut(FUNCTION_NAME);
}


void Foam::unsteadyTimeManipulation::storeTime()
{
    storedTime_.reset(new dimensionedScalar(time_));
    storedTimeIndex_.reset(new label(time_.timeIndex()));
    storedEndTime_.reset(new dimensionedScalar(time_.endTime()));
}


void Foam::unsteadyTimeManipulation::restoreTime()
{
    if (storedTime_ && storedTimeIndex_)
    {
        time_.setTime(storedTime_(), storedTimeIndex_());
        storedTime_.clear();
        storedTimeIndex_.clear();
    }
    else
    {
        FatalErrorInFunction
            << "Attempted to restore unset stored time and timeIndex"
            << exit(FatalError)
            << endl;
    }

    if (storedEndTime_)
    {
        time_.setEndTime(storedEndTime_());
        storedEndTime_.clear();
    }
    else
    {
        FatalErrorInFunction
            << "Attempted to restore unset stored endTime"
            << exit(FatalError)
            << endl;
    }
    printOut(FUNCTION_NAME);
}


void Foam::unsteadyTimeManipulation::writeEntries(Ostream& os) const
{
    os.beginBlock("localEntries");
    os.writeEntry("startTime",  startTime());
    os.writeEntry("startTimeIndex", startTimeIndex());
    os.writeEntry("primalEndTime", primalEndTime());
    os.writeEntry("endTimeIndex", endTimeIndex());
    os.writeEntry("span", span_);
    os.endBlock();

    os.beginBlock("TimeEntries");
    os.writeEntry("time", time_.value());
    os.writeEntry("timeIndex", time_.timeIndex());
    os.writeEntry("deltaT", time_.deltaT());
    os.writeEntry("endTime", time_.endTime());
    os.endBlock();
}


bool Foam::unsteadyTimeManipulation::movePoints()
{
    // Does nothing
    return true;
}


void Foam::unsteadyTimeManipulation::updateMesh(const mapPolyMesh&)
{
    // Does nothing
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const unsteadyTimeManipulation& timeManip
)
{
    timeManip.writeEntries(os);
    return os;
}


// ************************************************************************* //
