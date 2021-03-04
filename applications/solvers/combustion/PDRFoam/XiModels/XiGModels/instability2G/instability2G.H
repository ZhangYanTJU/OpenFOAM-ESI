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
    Foam::XiGModels::instability2G

Description
    Flame-surface instabilityG flame-wrinkling generation rate coefficient model
    used in \link XiModel.H \endlink.

SourceFiles
    instability2G.C

\*---------------------------------------------------------------------------*/

#ifndef instability2G_H
#define instability2G_H

#include "laminarFlameSpeed.H"
#include "SCOPELaminarFlameSpeed.H"
#include "XiGModel.H"
#include "ignitionSite.H"
#include "ignition.H"
#include "Time.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiGModels
{

/*---------------------------------------------------------------------------*\
                          Class instability2G Declaration
\*---------------------------------------------------------------------------*/

class instability2G
:
    public XiGModel
{
    // Private Data

        laminarFlameSpeedModels::SCOPE saModel_;

        // GIn (inituial instability G)divided by Su^2. Read from fuel file
        scalar CIn_;

        //- Default CIn if not in fuel file
        scalar defaultCIn_;

        // Determine how fast GIn fades out as turbulence starts
        scalar GInFade_;

        // Set GIn large so that XiEq determines Xi value.
        // Son increase byfactor:
        scalar GInMult_;

        //- The instability2G length-scale
        dimensionedScalar lambdaIn_;

        //- Xi generation rate model due to all other processes
        autoPtr<XiGModel> XiGModel_;


    // Private Member Functions

        //- No copy construct
        instability2G(const instability2G&) = delete;

        //- No copy assignment
        void operator=(const instability2G&) = delete;


public:

    //- Runtime type information
    TypeName("instability2G");


    // Constructors

        //- Construct from components
        instability2G
        (
            const dictionary& XiGProperties,
            const word& modelType,
            const psiuReactionThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su
        );


    //- Destructor
    virtual ~instability2G() = default;


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
