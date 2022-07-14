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


Class
    Foam::unsteadyTimeManipulation

Description
    Class with helper functions to manipulate time during unsteady optimisation



SourceFiles
    unsteadyTimeManipulation.C

\*---------------------------------------------------------------------------*/

#ifndef unsteadyTimeManipulation_H
#define unsteadyTimeManipulation_H

#include "MeshObject.H"
#include "fvMesh.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class unsteadyTimeManipulation;
Ostream& operator<< (Ostream&, const unsteadyTimeManipulation&);


/*---------------------------------------------------------------------------*\
                   Class unsteadyTimeManipulation Declaration
\*---------------------------------------------------------------------------*/

class unsteadyTimeManipulation
:
    public MeshObject<fvMesh, UpdateableMeshObject, unsteadyTimeManipulation>
{
protected:


    // Protected data

        //- Not constant reference to Time
        Time& time_;

        //- Start time of the optimisation loop
        //  Will be updated at the beginning of each new loop
        dimensionedScalar startTime_;

        //- Start time index of the optimisation loop
        label startTimeIndex_;

        //- Start time index of the optimisation loop
        label endTimeIndex_;

        //- Span of the optimisation loop
        scalar span_;


        // Temporarily stored time features

            autoPtr<dimensionedScalar> storedTime_;
            autoPtr<label> storedTimeIndex_;
            autoPtr<dimensionedScalar> storedEndTime_;


    // Protected Member Functions

        //- Print entries for debugging
        void printOut(const char* functionName);


private:

    // Private Member Functions

        //- No copy construct
        unsteadyTimeManipulation(const unsteadyTimeManipulation&) = delete;

        //- No copy assignment
        void operator=(const unsteadyTimeManipulation&) = delete;


public:

    //- Runtime type information
    TypeName("unsteadyTimeManipulation");


    // Constructors

        //- Construct from components
        unsteadyTimeManipulation(const fvMesh& mesh);


    //- Destructor
    virtual ~unsteadyTimeManipulation() = default;


    // Member Functions

        // Access

            //- Return start time of the optimisation loop
            inline const dimensionedScalar& startTime() const
            {
                return startTime_;
            }

            //- Return start time of the optimisation loop
            inline label startTimeIndex() const
            {
                return startTimeIndex_;
            }

            //- Return span of the optimisation loop
            inline label endTimeIndex() const
            {
                return endTimeIndex_;
            }

            //- Return span of the optimisation loop
            inline scalar span() const
            {
                return span_;
            }

            //- Return end time of primal
            inline scalar primalEndTime() const
            {
                return startTime_.value() + span_;
            }


        // Edit

            //- Set startTime
            inline void setStartTime(const scalar startTime)
            {
                startTime_.value() = startTime;
            }

            //- Set startTimeIndex
            inline void setStartTimeIndex(const label startTimeIndex)
            {
                startTimeIndex_ = startTimeIndex;
            }

            //- Set span
            inline void setSpan(const scalar span)
            {
                span_ = span;
            }

            //- Set endTimeIndex
            inline void setEndTimeIndex(const label endTimeIndex)
            {
                endTimeIndex_ = endTimeIndex;
            }


            //- Store current endTimeIndex to local copy
            inline void storeEndTimeIndex()
            {
                endTimeIndex_ = time_.timeIndex();
            }

            //- Move time to the beginning of the primal
            //  Optionally sets endTime equal to startTime + span
            void moveToPrimalStartTime(bool setEndTime = true);

            //- Move time to the end of the primal
            //  Optionally sets endTime equal to startTime
            void moveToAdjointStartTime
            (
                bool setEndTime = true,
                bool setTimeIndex = true
            );

            //- Increment time features at the beginning of each optimisation
            //- cycle
            void newOptimisationCycle();


        //- Store current time features
        void storeTime();

        //- Retrieve last stored time features
        void restoreTime();

        //- Write local entries and time entries to stream
        void writeEntries(Ostream& os) const;

        //- Dummy function required by MeshObject.
        virtual bool movePoints();

        //- Dummy function required by MeshObject.
        virtual void updateMesh(const mapPolyMesh&);


    // Ostream operators

        friend Ostream& operator<<
        (
            Ostream&,
            const unsteadyTimeManipulation&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
