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
    Foam::XiEqModel::BLMgMaXiEq

Description
    Model for XiEq based on Bradley, Lawes and Mansour (2011)
    Cobustion and Falme, 158, 123 correlation
    with a linear correction function to give a plausible profile for XiEq.
    See \link SCOPELaminarFlameSpeed.H \endlink for details on the SCOPE laminar
    flame speed model.

SourceFiles
    BLMgMaXiEq.C

\*---------------------------------------------------------------------------*/

#ifndef BLMgMaXiEq_H
#define BLMgMaXiEq_H

#include "XiEqModel.H"
#include "SCOPELaminarFlameSpeed.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiEqModels
{

/*---------------------------------------------------------------------------*\
                          Class BLMgMaXiEq Declaration
\*---------------------------------------------------------------------------*/

class BLMgMaXiEq
:
    public XiEqModel
{
    // Private data

        // Model constants

        scalar kaCoef_;
        scalar lowK0_;
        scalar lowKg_;
        scalar XiEqCoef_;
        scalar alphaCoefP_;
        scalar betaCoefP_;
        scalar alphaCoefN_;
        scalar betaCoefN_;
        scalar maLim_;
        scalar maLim1_;
        scalar quenchCoef_, quenchExp_, quenchM_;
        scalar quenchRate1_, quenchRate2_;
        scalar lCoef_;

        //- Minimum Su
        dimensionedScalar SuMin_;

        //- Schelkin effect Model constants
        scalar uPrimeCoef_;
        scalar nrExp_;

        //- Use sub-grid Schelkin effect
        bool subGridSchelkin_;

        //- The SCOPE laminar flame speed model used to obtain the
        //  Marstein number.  Note: the laminar flame speed need not be
        //  obtained form the same model.
        laminarFlameSpeedModels::SCOPE MaModel;


    // Private Member Functions

        //- No copy construct
        BLMgMaXiEq(const BLMgMaXiEq&) = delete;

        //- No copy assignment
        void operator=(const BLMgMaXiEq&) = delete;


public:

    //- Runtime type information
    TypeName("BLMgMaXiEq");


    // Constructors

        //- Construct from components
        BLMgMaXiEq
        (
            const dictionary& XiEqProperties,
            const word& modelType,
            const psiuReactionThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su
        );


    //- Destructor
    virtual ~BLMgMaXiEq();


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
