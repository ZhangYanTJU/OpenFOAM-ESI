/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::cyclicPeriodicAMILduInterface

Description
    An abstract base class for cyclic periodic AMI coupled interfaces

SourceFiles
    cyclicPeriodicAMILduInterface.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicPeriodicAMILduInterface_H
#define cyclicPeriodicAMILduInterface_H

#include "cyclicAMILduInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class cyclicPeriodicAMILduInterface Declaration
\*---------------------------------------------------------------------------*/

class cyclicPeriodicAMILduInterface
:
    public cyclicAMILduInterface
{
public:

    //- Runtime type information
    TypeName("cyclicPeriodicAMILduInterface");


    // Constructors

        //- Construct null
        cyclicPeriodicAMILduInterface() = default;


    //- Destructor
    virtual ~cyclicPeriodicAMILduInterface() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
