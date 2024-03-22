/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::fv::treeTurbulence

Group
    grpFvOptionsSources

Description
    Applies sources on \c k and either \c epsilon or \c omega
    to incorporate effects of trees on atmospheric boundary layer.

    Corrections applied to:
    \verbatim
      k         | Turbulent kinetic energy                  [m^2/s^2]
    \endverbatim

    Corrections applied to either of the below, if exist:
    \verbatim
      epsilon   | Turbulent kinetic energy dissipation rate [m^2/s^3]
      omega     | Specific dissipation rate                 [1/s]
    \endverbatim

    Required fields:
    \verbatim
      k         | Turbulent kinetic energy                  [m^2/s^2]
      Cd        | Canopy drag coefficient                   [-]
      LAD       | Leaf area density                         [m^2/m^3]
      epsilon   | Turbulent kinetic energy dissipation rate [m^2/s^3]
      omega     | Specific dissipation rate                 [1/s]
    \endverbatim

    References:
    \verbatim
        Governing equations (tag:BSG):
            Brozovsky, J., Simonsen, A., & Gaitani, N. (2021).
            Validation of a CFD model for the evaluation of urban microclimate
            at high latitudes: A case study in Trondheim, Norway.
            Building and Environment, 205, 108175.
            DOI:10.1016/j.buildenv.2021.108175
    \endverbatim

Usage
    Example by using \c constant/fvOptions:
    \verbatim
    treeTurbulence1
    {
        // Mandatory entries
        type         treeTurbulence;

        // Optional entries
        Cd           <word>;
        LAD          <word>;
        betaP        <scalar>;
        betaD        <scalar>;
        Ceps1        <scalar>;
        Ceps2        <scalar>;
        betaStar     <scalar>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                         | Type | Reqd | Deflt
      type     | Type name: treeTurbulence           | word | yes  | -
      Cd       | Name of operand canopy drag coefficient field | word | no | Cd
      LAD      | Name of operand leaf area density field | word | no | LAD
      betaP    | Model coefficient                 | scalar | no   | 1.0
      betaD    | Model coefficient                 | scalar | no   | 4.0
      Ceps1    | Model coefficient                 | scalar | no   | 0.9
      Ceps2    | Model coefficient                 | scalar | no   | 0.9
      betaStar | Model coefficient                 | scalar | no   | 0.09
    \endtable

    The inherited entries are elaborated in:
      - \link fvOption.H \endlink
      - \link cellSetOption.H \endlink

SourceFiles
    treeTurbulence.C
    treeTurbulenceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fv_treeTurbulence_H
#define fv_treeTurbulence_H

#include "cellSetOption.H"
#include "turbulentTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class treeTurbulence Declaration
\*---------------------------------------------------------------------------*/

class treeTurbulence
:
    public fv::cellSetOption
{
    // Private Data

        //- Internal flag to determine the working field is epsilon or omega
        bool isEpsilon_;

        // Model coefficients

            scalar betaP_;
            scalar betaD_;
            scalar Ceps1_;
            scalar Ceps2_;
            scalar betaStar_;

        //- Name of operand canopy drag coefficient field
        word CdName_;

        //- Name of operand leaf area density field
        word LADname_;


     // Private Member Functions

        //- Return requested field from the object registry
        //- or read+register the field to the object registry
        volScalarField& getOrReadField(const word& fieldName) const;

        //- Apply sources to turbulent kinetic energy
        template<class AlphaFieldType, class RhoFieldType>
        void kSource
        (
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            fvScalarMatrix& eqn,
            const label fieldi
        ) const;

        //- Apply sources to turbulent kinetic energy dissipation rate
        template<class AlphaFieldType, class RhoFieldType>
        void epsilonSource
        (
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            fvScalarMatrix& eqn,
            const label fieldi
        ) const;

        //- Apply sources to specific dissipation rate
        template<class AlphaFieldType, class RhoFieldType>
        void omegaSource
        (
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            fvScalarMatrix& eqn,
            const label fieldi
        ) const;


public:

    //- Runtime type information
    TypeName("treeTurbulence");


    // Constructors

        //- Construct from explicit source name and mesh
        treeTurbulence
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        treeTurbulence(const treeTurbulence&) = delete;

        //- No copy assignment
        void operator=(const treeTurbulence&) = delete;


    // Member Functions

        //- Add explicit contribution to epsilon or omega equation
        //- for incompressible flow computations
         virtual void addSup
        (
            fvScalarMatrix& eqn,
            const label fieldi
        );

        //- Add explicit contribution to epsilon or omega equation
        //- for compressible flow computations
        virtual void addSup
        (
            const volScalarField& rho,
            fvScalarMatrix& eqn,
            const label fieldi
        );

        //- Add explicit contribution to epsilon or omega equation
        //- for multiphase flow computations
        virtual void addSup
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            fvScalarMatrix& eqn,
            const label fieldi
        );

        //- Read source dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "treeTurbulenceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
