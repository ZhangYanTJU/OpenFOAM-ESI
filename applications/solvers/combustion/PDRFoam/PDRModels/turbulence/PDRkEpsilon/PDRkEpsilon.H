/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::compressible::RASModels::PDRkEpsilon

Description
    Standard k-epsilon turbulence model with additional source terms
    corresponding to PDR basic drag model (\link basic.H \endlink)

SourceFiles
    PDRkEpsilon.C

\*---------------------------------------------------------------------------*/

#ifndef compressiblePDRkEpsilon_H
#define compressiblePDRkEpsilon_H

#include "kEpsilon.H"
#include "ignition.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class PDRkEpsilon Declaration
\*---------------------------------------------------------------------------*/

class PDRkEpsilon
:
    public Foam::RASModels::kEpsilon
    <
        EddyDiffusivity
        <
            compressible::turbulenceModel
        >
    >
{
    // Private Data

        // Model coefficients
        scalar C5_;
        scalar C6_;
        scalar maxLOverCellW_;
        scalar lCoef_;
        dimensionedScalar noTurbUntil_;
        dimensionedScalar LOverLobs_;
        dimensionedScalar LOverMobs_;
        scalar bMin_;
        bool flameFilteredG_;


public:

    //- Runtime type information
    TypeName("PDRkEpsilon");


    // Constructors

        //- Construct from components
        PDRkEpsilon
        (
            const geometricOneField& alpha,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const fluidThermo& thermophysicalModel,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~PDRkEpsilon();


    // Member Functions

        //- Solve the turbulence equations and correct the turbulence viscosity
        void correct();

        //- Read turbulenceProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
