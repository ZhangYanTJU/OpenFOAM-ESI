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
    Foam::XiEqModels::normBasicSubGrid

Description


SourceFiles
    normBasicSubGrid.C

\*---------------------------------------------------------------------------*/

#ifndef normBasicSubGrid_H
#define normBasicSubGrid_H

#include "XiEqModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiEqModels
{

/*---------------------------------------------------------------------------*\
                          Class normBasicSubGrid Declaration
\*---------------------------------------------------------------------------*/

class normBasicSubGrid
:
    public XiEqModel
{
    // Private Data

        // Constants in the equilibrium Xp equation
        scalar Cxpe1_;
        scalar Cxpe2_;
        scalar Cxpe3_;
        scalar Cxpe4_;


    // Private Member Functions

        //- No copy construct
        normBasicSubGrid(const normBasicSubGrid&) = delete;

        //- No copy assignment
        void operator=(const normBasicSubGrid&) = delete;


public:

    //- Runtime type information
    TypeName("normBasicSubGrid");


    // Constructors

        //- Construct from components
        normBasicSubGrid
        (
            const dictionary& XiEqProperties,
            const word& modelType,
            const psiuReactionThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su
        );


    //- Destructor
    virtual ~normBasicSubGrid();


    // Member Functions

        //- Return the flame-wrinkling XiEq
        virtual tmp<volScalarField> XiEq() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& XiEqProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiEqModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
