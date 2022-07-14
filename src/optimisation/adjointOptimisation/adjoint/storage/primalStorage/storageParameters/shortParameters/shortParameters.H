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
    shortParameters

Description
    Reads from dictionary and stores all parameters related to the storage
    of the primal fields

SourceFiles
    shortParameters.C

\*---------------------------------------------------------------------------*/

#ifndef shortParameters_H
#define shortParameters_H

#include "storageParameters.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class shortParameters Declaration
\*---------------------------------------------------------------------------*/

class shortParameters
:
    public storageParameters
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        shortParameters(const shortParameters&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const shortParameters&) = delete;


public:

    TypeName("shortParameters");

    // Constructors

        //- Construct from components
        shortParameters
        (
            const fvMesh& mesh,
            const dictionary& storageDict,
            const wordList& allocatedFieldNames
        );


    //- Destructor
    virtual ~shortParameters() = default;


    // Member functions

        virtual void initialize() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
