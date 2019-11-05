/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |   Copyright 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
   acousticFoam 

Group

Description
    Acoustic solver solving the acoustic pressure wave equation.
    
    \f[
            \ddt2{pa} - sqr(c) \laplacian{pa} = 0
    \f]
    
    Where:
    \vartable
        c       | Sound speed
        pa      | Acoustic pressure
    \endvartable

SourceFiles
    acousticFoam.C

company
    Volkswagen Group Research 
    by
    Alexander Kabat vel Job
    email:alexander.kabat.vel.job@volkswagen.de 


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Acoustic solver solving the acoustic pressure wave equation."
    );
    
    #include "postProcess.H"
    
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "readTransportProperties.H"
    #include "createFields.H"

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {   
        ++runTime;
        
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "paEqn.H"

        runTime.write();
        
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;
    
    return(0);
}


// ************************************************************************* //
