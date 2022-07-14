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
    Foam::PIMPLEControlOpt

Description
    PIMPLE control class for single runs (i.e. not optimisation).
    Time acts as in simpleFoam, with all solver control read through
    optimisationDict

\*---------------------------------------------------------------------------*/

#ifndef PIMPLEControlOpt_H
#define PIMPLEControlOpt_H

#include "PIMPLEControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class PIMPLEControlOpt Declaration
\*---------------------------------------------------------------------------*/

class PIMPLEControlOpt
:
    public PIMPLEControl
{
private:

        //- No copy construct
        PIMPLEControlOpt(const PIMPLEControlOpt&) = delete;

        //- No copy assignment
        void operator=(const PIMPLEControlOpt&) = delete;


public:


    // Static Data Members

        //- Run-time type information
        TypeName("unsteadyOptimisation");


    // Constructors

        //- Construct from mesh
        PIMPLEControlOpt
        (
            fvMesh& mesh,
            const word& managerType,
            const solver& solver
        );


    //- Destructor
    virtual ~PIMPLEControlOpt() = default;


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
