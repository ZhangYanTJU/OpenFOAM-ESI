/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::fv::buoyancyTurbSource

Group
    grpFvOptionsSources

Description
    Applies sources on \c k and/or \c omega to incorporate
    effects of buoyancy in incompressible and compressible flows.

    For incompressible cases:

    Corrections applied to:
    \verbatim
      k         | Turbulent kinetic energy                   [m2/s2]
      omega     | Specific dissipation rate                  [1/s]
    \endverbatim

    Required fields:
    \verbatim
      k         | Turbulent kinetic energy                   [m2/s2]
      omega     | Specific dissipation rate                  [1/s]
      alphat    | Kinematic turbulent thermal conductivity   [m2/s]
      T         | Temperature                                [K]
    \endverbatim


    For compressible cases:

    Corrections applied to:
    \verbatim
      k         | Turbulent kinetic energy                   [m2/s2]
    \endverbatim

    Required fields:
    \verbatim
      k         | Turbulent kinetic energy                   [m2/s2]
      rho       | Fluid density                              [kg/m3]
    \endverbatim

    References:
    \verbatim
        Buoyancy modifier for compressible flow cases (tag:DTR):
            Devolder, B., Troch, P., & Rauwoens, P. (2018).
            Performance of a buoyancy-modified k-ω and k-ω
            SST turbulence model for simulating wave breaking
            under regular waves using OpenFOAM.
            Coastal Engineering, 138, 49-65.
            DOI:10.1016/j.coastaleng.2018.04.011
    \endverbatim

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    buoyancyTurbSource1
    {
        // Mandatory entries (unmodifiable)
        type                 buoyancyTurbSource;

        buoyancyTurbSourceCoeffs
        {
            // Mandatory (inherited) entries (unmodifiable)
            selectionMode    all;

            // Optional entries (unmodifiable)
            rho             rho;
            beta            3.3e-03;
        }

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type   | Reqd  | Dflt
      type      | Type name: buoyancyTurbSource     | word   | yes   | -
      rho       | Name of density field             | word   | no    | rho
      beta      | Thermal expansion coefficient for incompressible cases <!--
                                             -->    | scalar | no    | 3.3e-03
    \endtable

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink
     - \link cellSetOption.H \endlink

Note
  - Implementation for incompressible cases follows
    a heuristic approach of an industrial partner.
  - Personal communication with Dr. Devolder for compressible cases (28 Oct 20):

    "Eq. 21 of my publication is the buoyancy source term. This term
    appears when the transport equation for the TKE is fundamentally
    derived from the Favre-averaged (density weighted) low Mach number
    equations, as reported in:

        Van Maele, K., & Merci, B. (2006).
        Application of two buoyancy-modified k–ε
        turbulence models to different types of buoyant plumes.
        Fire Safety Journal, 41(2), 122-138.
        DOI:10.1016/j.firesaf.2005.11.003

        Van Maele, K., & Merci, B. (2006).
        Importance of buoyancy and chemistry modelling in
        steady RANS simulations of well-ventilated tunnel fires.
        Turkish Jour. of Engineering and Environmental Sciences, 30(3), 145-155.

    Therefore, this term is generally applicable for
    applications where density gradients exist such as
    compressible heat transfer but also wave modelling
    due to the VoF method."

    "Due to VoF, the interface between water and air is not
    a jump but it is rather smooth. This means that there is
    a finite density gradient and therefore I think that the
    density must be included inside the differential operators."

SourceFiles
    buoyancyTurbSource.C
    buoyancyTurbSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fv_buoyancyTurbSource_H
#define fv_buoyancyTurbSource_H

#include "cellSetOption.H"
#include "turbulentTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class buoyancyTurbSource Declaration
\*---------------------------------------------------------------------------*/

class buoyancyTurbSource
:
    public cellSetOption
{
    // Private Data

        //- Name of density field
        const word rhoName_;

        //- Thermal expansion coefficient for incompressible cases [1/K]
        const dimensionedScalar beta_;


        // Fields

            //- Buoyancy production term for incompressible cases [m2/s3]
            volScalarField::Internal B_;


     // Private Member Functions

        //- Compute the buoyancy production term for incompressible cases
        void calcB();

        //- Apply buoyancyTurbSource to omega field for incompressible cases
        void buoyancyTurbSourceOmega(fvMatrix<scalar>& eqn) const;

        //- Apply buoyancyTurbSource to k field for incompressible cases
        void buoyancyTurbSourceK(fvMatrix<scalar>& eqn) const;

        //- Apply buoyantTurbSource to k field for compressible cases
        template<class AlphaFieldType, class RhoFieldType>
        void buoyancyTurbSourceK
        (
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        ) const;

public:

    //- Runtime type information
    TypeName("buoyancyTurbSource");


    // Constructors

        //- No default construct
        buoyancyTurbSource() = delete;

        //- Construct from explicit source name and mesh
        buoyancyTurbSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        buoyancyTurbSource(const buoyancyTurbSource&) = delete;

        //- No copy assignment
        void operator=(const buoyancyTurbSource&) = delete;


    // Member Functions

        //- Add explicit contribution to omega equation
        //- for incompressible flow computations
        virtual void addSup
        (
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Add explicit contribution to omega equation
        //- for compressible flow computations
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Add explicit contribution to omega equation
        //- for multiphase flow computations
        virtual void addSup
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Read source dictionary (effectively no-op)
        virtual bool read(const dictionary& dict)
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "buoyancyTurbSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
