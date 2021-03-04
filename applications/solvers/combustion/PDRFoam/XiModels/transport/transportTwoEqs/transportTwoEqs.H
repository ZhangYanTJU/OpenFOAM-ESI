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
    Foam::XiModels::transportTwoEqs

Description

SourceFiles
    transportTwoEqs.C

\*---------------------------------------------------------------------------*/

#ifndef transportTwoEqs_H
#define transportTwoEqs_H

#include "XiModel.H"
#include "XiEqModel.H"
#include "XiGModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiModels
{

/*---------------------------------------------------------------------------*\
                          Class transportTwoEqs Declaration
\*---------------------------------------------------------------------------*/

class transportTwoEqs
:
    public XiModel
{
    // Private data

        scalar XiShapeCoef_;
        scalar CpfiDot_;
        scalar CpfiCross_;
        scalar GEtaExp_;

        //- Maximum Lobs/CellWidth
        scalar LOverCw_;

        //- Equilibrium for Xi (turbulence)
        autoPtr<XiEqModel> XiEqModel_;

        //- Generation for Xi (turbulence)
        autoPtr<XiGModel> XiGModel_;

        //- Equilibrium for Xp (obstacles)
        autoPtr<XiEqModel> XpEqModel_;

        //- Generation for Xp (obstacles)
        autoPtr<XiGModel> XpGModel_;

        //- Dissipation length scale for subgrid obstacles
        volScalarField Ep_;


    // Private Member Functions

        //- No copy construct
        transportTwoEqs(const transportTwoEqs&) = delete;

        //- No copy assignment
        void operator=(const transportTwoEqs&) = delete;


public:

    //- Runtime type information
    TypeName("transportTwoEqs");


    // Constructors

        //- Construct from components
        transportTwoEqs
        (
            const dictionary& XiProperties,
            const psiuReactionThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su,
            const volScalarField& rho,
            const volScalarField& b,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~transportTwoEqs();


    // Member Functions

        //- Return the flame diffusivity
        virtual tmp<volScalarField> Db() const;

        //- Correct the flame-wrinkling Xi
        virtual void correct()
        {
            NotImplemented;
        }

        //- Correct the flame-wrinkling Xi using the given convection scheme
        virtual void correct(const fv::convectionScheme<scalar>& mvConvection);

        //- Update properties from given dictionary
        virtual bool read(const dictionary& XiProperties);

        //- Write fields of the XiEq model
        virtual void writeFields()
        {
            XiEqModel_().writeFields();
            XpEqModel_().writeFields();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
