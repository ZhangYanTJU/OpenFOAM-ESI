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
    Foam::XiEqModels::instability2XiEq

Description


SourceFiles
    instability2XiEq.C

\*---------------------------------------------------------------------------*/

#ifndef instability2XiEq_H
#define instability2XiEq_H

#include "laminarFlameSpeed.H"
#include "SCOPELaminarFlameSpeed.H"
#include "ignitionSite.H"
#include "ignition.H"
#include "Time.H"
#include "fvMesh.H"
#include "XiEqModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiEqModels
{

/*---------------------------------------------------------------------------*\
                          Class instability2 Declaration
\*---------------------------------------------------------------------------*/

class instability2XiEq
:
    public XiEqModel
{
    // Private data

        // -Laminar burning speed
        laminarFlameSpeedModels::SCOPE saModel_;

        //- GIn (initial instability G)divided by Su^2. Read from fuel file
        scalar CIn_;

        //- Default CIn if not in fuel file
        scalar defaultCIn_;

        //- Determines how fast XiEqIn fades out as turbulence comes in
        scalar XiEqInFade_;

        //- Equilibrium Xi model due to all other effects
        autoPtr<XiEqModel> XiEqModel_;


    // Private Member Functions

        //- Disallow copy construct
        instability2XiEq(const instability2XiEq&);

        //- Disallow default bitwise assignment
        void operator=(const instability2XiEq&);


public:

    //- Runtime type information
    TypeName("instability2XiEq");


    // Constructors

        //- Construct from components
        instability2XiEq
        (
            const dictionary& XiEqProperties,
            const word& modelType,
            const psiuReactionThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su
        );


    //- Destructor
    virtual ~instability2XiEq();


    // Member Functions

        //- Return the flame-wrinking XiEq
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
