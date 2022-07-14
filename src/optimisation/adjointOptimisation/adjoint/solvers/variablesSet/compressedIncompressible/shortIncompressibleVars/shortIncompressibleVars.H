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
    Foam::shortIncompressibleVars

Description
    Compresses the current time-instance known by the primal solvers.
    Does not include oldTimes

\*---------------------------------------------------------------------------*/

#ifndef shortIncompressibleVars_H
#define shortIncompressibleVars_H

#include "compressedIncompressibleVars.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class shortIncompressibleVars Declaration
\*---------------------------------------------------------------------------*/

class shortIncompressibleVars
:
    public compressedIncompressibleVars
{
public:


    // Static Data Members

        //- Run-time type information
        TypeName("shortIncompressibleVars");


    // Constructors

        shortIncompressibleVars
        (
            incompressibleVars& vs,
            storageParameters& storageParams
        );


    //- Destructor
    virtual ~shortIncompressibleVars() = default;


    // Member Functions

        //- Compress the incompressibleVars
        void compress() override;

        //- Function to return the compression metrics (initial size,
        //- compressed size if compression is applied, etc)
        void decompress
        (
            incompressibleVars& vars
        ) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
