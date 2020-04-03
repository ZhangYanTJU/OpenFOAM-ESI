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

Class
    Foam::functionObjects::setTimeStepFaRegionsFunctionObject

Group
    grpUtilitiesFunctionObjects

Description
    This function object control the time step for classes of the type
    regionFaModel. It reads regionFaMaxCo entro from controlDict and
    evaluate the timne step based on the finite area Courant Number.

    Can only be used with solvers with adjustTimeStep control (e.g.
    pimpleFoam).  It makes no attempt to co-operate with other time step
    'controllers', e.g. maxCo, other functionObjects. Supports 'enabled'
    flag but none of the other options 'timeStart', 'timeEnd', 'writeControl'
    etc.

Usage
    Example of function object specification to manipulate the time step:
    \verbatim
    setTimeStep1
    {
        type        setTimeStepFaRegion;

        ...
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | setTimeStepFaRegion     | yes         |
        enabled      | On/off switch           | no          | yes
        timeStart    | Start time              | no          | 0
        timeEnd      | End time                | no          | GREAT
    \endtable

SourceFiles
    setTimeStepFaRegionsFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_setTimeStepFaRegionsFunctionObject_H
#define functionObjects_setTimeStepFaRegionsFunctionObject_H

#include "timeFunctionObject.H"
#include "regionFaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam::regionModels;

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                Class setTimeStepFaRegionsFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class setTimeStepFaRegionsFunctionObject
:
    public functionObjects::timeFunctionObject
{
    // Private Member Data

        //- Time step function/table
        //autoPtr<Function1<scalar>> timeStepPtr_;


    // Private Member Functions

        //- No copy construct
        setTimeStepFaRegionsFunctionObject
        (
            const setTimeStepFaRegionsFunctionObject&
        ) = delete;

        //- No copy assignment
        void operator=(const setTimeStepFaRegionsFunctionObject&) = delete;

        //- Return minimum deltaT from fa regions
        scalar regionDeltaT() const;


public:

    //- Runtime type information
    TypeName("setTimeStepFaRegion");


    // Constructors

        //- Construct from components
        setTimeStepFaRegionsFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    // Destructor
    virtual ~setTimeStepFaRegionsFunctionObject() = default;


    // Member Functions

        //- Called at the end of Time::adjustDeltaT() if adjustTime is true
        virtual bool adjustTimeStep();

        //- Read and set the function object if its data have changed
        virtual bool read(const dictionary& dict);

        //- Execute does nothing
        virtual bool execute();

        //- Write does nothing
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
