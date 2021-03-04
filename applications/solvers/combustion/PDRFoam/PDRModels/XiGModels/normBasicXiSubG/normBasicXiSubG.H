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
    Foam::XiGModel::normBasicSubGrid

Description


SourceFiles
    normBasicSubGrid.C

\*---------------------------------------------------------------------------*/

#ifndef normBasicSubGrid_H
#define normBasicSubGrid_H

#include "XiGModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiGModels
{

/*---------------------------------------------------------------------------*\
                          Class normBasicSubGrid Declaration
\*---------------------------------------------------------------------------*/

class normBasicSubGrid
:
    public XiGModel
{
    // Private Data

        //- Sub-grid generation rate coefficient
        scalar k1_;

        //- Sub-grid generation rate coefficient - * sqrt(b)
        scalar kb1_;

        //- Sub-grid generation rate coefficient - * b
        scalar kbe_;

        //- Sub-grid generation rate upper limit coefficient - * b
        scalar kbx_;

        //- Sub-grid generation rate coefficient
        scalar k2_;

        //- Maximum Lobs/CellWidth
        scalar LOverCw_;


    // Private Member Functions

        //- No copy construct
        normBasicSubGrid(const normBasicSubGrid&) = delete;

        //- No copy assignment
        void operator=(const normBasicSubGrid&) = delete;


public:

    //- Runtime type information
    TypeName("normBasicSubGridG");


    // Constructors

        //- Construct from components
        normBasicSubGrid
        (
            const dictionary& XiGProperties,
            const word& modelType,
            const psiuReactionThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su
        );


    //- Destructor
    virtual ~normBasicSubGrid() = default;


    // Member Functions

        //- Return the flame-wrinkling generation rate
        virtual tmp<volScalarField> G() const;

        //- Return the flame diffusivity
        virtual tmp<volScalarField> Db() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& XiGProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiGModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
