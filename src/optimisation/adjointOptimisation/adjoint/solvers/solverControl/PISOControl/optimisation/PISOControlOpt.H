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
    Foam::PISOControlOpt

Description
    PISO control class for single runs (i.e. not optimisation).
    Time acts as in simpleFoam, with all solver control read through
    optimisationDict

\*---------------------------------------------------------------------------*/

#ifndef PISOControlOpt_H
#define PISOControlOpt_H

#include "PISOControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class PISOControlOpt Declaration
\*---------------------------------------------------------------------------*/

class PISOControlOpt
:
    public PISOControl
{
private:

        //- No copy construct
        PISOControlOpt(const PISOControlOpt&) = delete;

        //- No copy assignment
        void operator=(const PISOControlOpt&) = delete;


public:


    // Static Data Members

        //- Run-time type information
        TypeName("unsteadyOptimisation");


    // Constructors

        //- Construct from mesh
        PISOControlOpt
        (
            fvMesh& mesh,
            const word& managerType,
            const solver& solver
        );


    //- Destructor
    virtual ~PISOControlOpt() = default;


    // Member Functions

        // Solution control

            //- Whether to call time.write() or not
            virtual void writeFields();


        // Evolution

            //- Loop
            virtual bool loop();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
