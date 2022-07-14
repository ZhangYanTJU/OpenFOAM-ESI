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

Class
    checkPoint

Description

SourceFiles
    checkPoint.C

\*---------------------------------------------------------------------------*/

#ifndef checkPoint_H
#define checkPoint_H

#include "fullIncompressibleVars.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class checkPoint Declaration
\*---------------------------------------------------------------------------*/

class checkPoint
{
protected:

    // Protected data

        //- Variables held by this checkPoint
        autoPtr<fullIncompressibleVars> vars_;

        storageParameters& storageParams_;

        const fvMesh& mesh_;

        scalarList storageMetrics_;

        //- deltaT corresponding to the time-step this checkPoint holds
        scalar deltaT_;

        //- Whether this checkPoint is active
        Switch active_;

        //- Whether this checkPoint is a placeHolder
        Switch placeHolder_;

        //- schedule with size equal to nTimeSteps containing the number of
        //- time-solution that must be stored at the checkPoint:
        //  -1: inactive checkPoint     0: placeholder checkPoint;
        //   1: single checkPoint       2: double checkPoint
        label a_;

        label checkPointLevel_;

        label checkPointTimeIndex_;

        scalar checkPointTime_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        checkPoint(const checkPoint&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const checkPoint&) = delete;

public:

    // Static Data Members

        //- Run-time type information
        TypeName("checkPoint");

    // Constructors

        checkPoint(storageParameters& storageParams);


    //- Destructor
    virtual ~checkPoint() = default;


    // Access functions

        inline scalar deltaT() const
        {
            return deltaT_;
        }

        inline scalar& deltaT()
        {
            return deltaT_;
        }

        inline Switch active() const
        {
            return active_;
        }

        inline Switch& active()
        {
            return active_;
        }

        inline const Switch& placeHolder() const
        {
            return placeHolder_;
        }

        inline Switch& placeHolder()
        {
            return placeHolder_;
        }

        inline label a() const
        {
            return a_;
        }

        inline label& a()
        {
            return a_;
        }

        inline label level() const
        {
            return checkPointLevel_;
        }

        inline label& level()
        {
            return checkPointLevel_;
        }

        inline label timeIndex() const
        {
            return checkPointTimeIndex_;
        }

        inline label& timeIndex()
        {
            return checkPointTimeIndex_;
        }

        inline scalar timeValue() const
        {
            return checkPointTime_;
        }

        inline scalar& timeValue()
        {
            return checkPointTime_;
        }


    // Member functions

        //- Function to decompress the incompressible vars and update the
        //- solver fields with the decompressed ones
        void retrieve();

        void setPlaceHolder(Switch newLevel);

        //- store function must be called after setPlaceHolder.
        void store(incompressibleVars& vs, const label a = 2);

        //- Clear the chekpoint memory and reset its properties
        void empty();

        //- Clear the storage of the checkpoint and make it a placeHolder
        //- checkpoint
        void resetToPlaceHolder();

        const scalarList& storageMetrics() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
