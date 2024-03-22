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
    Foam::evapotranspirationHeatTransferModels::grass

Description
    Applies sources on temperature (\c T - incompressible) or energy
    (\c h/e - compressible) equation to incorporate evapotranspiration
    heat-transfer effects from the specified grass canopy.
    The model is based on Pemnan-Monteith Equation model.

    Sources applied to either of the below, if exist:
    \verbatim
      T         | Temperature                                [K]
      e         | Internal energy                            [m^2/s^2]
      h         | Enthalphy                                  [m^2/s^2]
    \endverbatim

    Required fields:
    \verbatim
      T         | Temperature                                [K]
      e         | Internal energy                            [m^2/s^2]
      h         | Enthalphy                                  [m^2/s^2]
      LAD       | Leaf area density                          [m^2/m^3]
    \endverbatim

    References:
    \verbatim
        Governing equations (tag:BSG):
            Brozovsky, J., Simonsen, A., & Gaitani, N. (2021).
            Validation of a CFD model for the evaluation of urban
            microclimate at high latitudes: A case study in Trondheim, Norway.
            Building and Environment, 205, 108175.
            DOI:10.1016/j.buildenv.2021.108175

        Governing equations (tag:APR):
            Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
            Crop evapotranspiration-Guidelines for computing crop
            water requirements-FAO Irrigation and drainage paper 56.
            Fao, Rome, 300(9), D05109.
    \endverbatim

Usage
    Example by using \c constant/fvOptions:
    \verbatim
    evapotranspirationHeatTransfer1
    {
        // Inherited entries
        ...

        // Mandatory entries
        model               grass;
        Tref                <Function1<scalar>>;
        RHref               <Function1<scalar>>;

        // Conditional entries

            // when soilHeatFluxMethod == boundary
            patches             <wordRes>;

        // Optional entries
        soilHeatFluxMethod  <word>;
        Csoil               <scalar>;
        rho                 <scalar>;
        Cp                  <scalar>;
        epsilon             <scalar>;
        h                   <scalar>;
        kappa               <scalar>;
        uRef                <scalar>;
        zRefU               <scalar>;
        zRefH               <scalar>;
        zom                 <scalar>;
        zoh                 <scalar>;
        d                   <scalar>;
        ri                  <scalar>;
        pAtm                <scalar>;
        T                   <word>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                         | Type | Reqd | Deflt
      model    | Model name: grass                   | word | yes  | -
      Tref     | Reference weather station air temperature [Celsius] <!--
               -->                     | Function1\<scalar\> | yes | -
      RHref    | Reference weather station relative humidity [%]     <!--
               -->                     | Function1\<scalar\> | yes | -
      patches  | Names of ground patches           | wordRes | yes | -
      soilHeatFluxMethod | Method to calculate soil heat flux - see below <!--
               -->                   | word | no | -
      Csoil    | Proportionality constant of soil heat-flux wrt <!--
               --> solar radiation                   | scalar | no | 0.1
      rho      | Mean air density at constant pressure [kg/m^3]      <!--
               -->                                   | scalar | no | 1.225
      Cp       | Specific heat at constant pressure [J/kg/C]         <!--
               -->                                   | scalar | no | 1013.0
      epsilon  | Molecular-weight ratio of water vapour/dry air [-]  <!--
               -->                                  | scalar | no  | 0.622
      h        | Height of grass layer [m]          | scalar | no  | 0.1
      kappa    | Von Karman constant [-]            | scalar | no  | 0.41
      uRef     | Reference velocity magnitude [m/s] | scalar | no  | 2.0
      zRefU    | Height of wind speed measurements [m] | scalar | no | 10.0
      zRefH    | Height of humidity measurements [m]   | scalar | no | 10.0
      zom      | Roughness length governing momentum transfer coefficient <!--
               -->                                  [-] | scalar | no | 0.123
      zoh      | Roughness length governing transfer of heat and vapour <!--
               --> coefficient [-]                      | scalar | no | 0.0123
      d        | Zero plane displacement height coefficient [-]         <!--
               -->                                      | scalar | no | 0.666
      ri       | Bulk stomata resistance of leaf [s/m] | scalar | no | 100.0
      pAtm     | Atmospheric pressure [kPa]            | scalar | no | 101.325
      T        | Name of temperature field             | word   | no | T
    \endtable

    Options for the \c soilHeatFluxMethod entry:
    \verbatim
      proportionalToSolarRadiation | Estimate from solar load
      boundary                     | Obtain soil heat flux from boundary
    \endverbatim

    The inherited entries are elaborated in:
      - \link evapotranspirationHeatTransfer.H \endlink
      - \link evapotranspirationHeatTransferModel.H \endlink

SourceFiles
    grass.C
    grassTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_evapotranspirationHeatTransferModels_grass_H
#define Foam_evapotranspirationHeatTransferModels_grass_H

#include "evapotranspirationHeatTransferModel.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace evapotranspirationHeatTransferModels
{

/*---------------------------------------------------------------------------*\
                          Class grass Declaration
\*---------------------------------------------------------------------------*/

class grass
:
    public evapotranspirationHeatTransferModel
{
    // Private Enumerations

        //- Options for the soil heat flux
        enum soilHeatFluxType : char
        {
            PROPORTIONAL_TO_SOLAR_RADIATION = 0, //!< "Estimate from solar load"
            BOUNDARY,            //!< "Obtain soil heat flux from boundary"
        };

        //- Names for soilHeatFluxType
        static const Enum<soilHeatFluxType> soilHeatFluxTypeNames;


    // Private Data

        //- Soil heat-flux calculation method
        enum soilHeatFluxType soilHeatFluxMethod_;

        //- Reference weather station air temperature [Celsius]
        autoPtr<Function1<scalar>> Tptr_;

        //- Reference weather station relative humidity [%]
        autoPtr<Function1<scalar>> RHptr_;

        //- Area coverage of grass [m^2]
        mutable scalar area_;

        //- Proportionality constant of soil heat-flux wrt solar radiation [-]
        scalar Csoil_;

        //- Mean air density at constant pressure [kg/m^3]
        scalar rho_;

        //- Specific heat at constant pressure [J/kg/C]
        scalar Cp_;

        //- Molecular-weight ratio of water vapour/dry air [-]
        scalar epsilon_;

        //- Height of grass layer [m]
        scalar h_;

        //- Von Karman constant [-]
        scalar kappa_;

        //- Reference velocity magnitude [m/s]
        scalar uRef_;

        //- Height of wind speed measurements [m]
        scalar zRefU_;

        //- Height of humidity measurements [m]
        scalar zRefH_;

        //- Roughness length governing momentum transfer coefficient [-]
        scalar zom_;

        //- Roughness length governing transfer of heat and vapour coeff [-]
        scalar zoh_;

        //- Zero plane displacement height coefficient [-]
        scalar d_;

        //- Bulk stomata resistance of leaf [s/m]
        scalar ri_;

        //- Atmospheric pressure [kPa]
        scalar pAtm_;

        //- Cached reference weather station air temperature [Celsius]
        mutable scalar Tref_;

        //- Name of temperature field
        word TName_;

        //- List of patches to calculate boundary heat flux
        labelHashSet patchSet_;


    // Private Member Functions

        //- Return evapotranspiration heat flux through grass layer [W/m^2]
        tmp<scalarField> E(const labelList& cells) const;

        //- Return slope of the relation between
        //- vapour pressure-temperature [Pa/K]
        scalar Delta() const;

        //- Return atmospheric vapour pressure deficit [Pa]
        scalar D() const;

        //- Return saturated vapour pressure over water [Pa]
        scalar pSat() const;

        //- Return psychrometric constant [kPa/K]
        scalar gamma() const;

        //- Return specific latent heat of vaporisation [MJ/kg]
        scalar lambda() const;

        //- Return bulk aerodynamic resistance [Pa]
        scalar ra() const;

        //- Return bulk surface resistance [Pa]
        scalar rs() const;

        //- Return area coverage of grass [m^2]
        scalar S(const labelList& cells) const;

        //- Return area-averaged heat flux through ground [W/m^2]
        scalar G(const labelList& cells) const;

        //- Return heat-flux boundary fields
        tmp<FieldField<Field, scalar>> qBf() const;


public:

    //- Runtime type information
    TypeName("grass");


    // Constructors

        //- Construct from components
        grass
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        grass(const grass&) = delete;

        //- No copy assignment
        void operator=(const grass&) = delete;


    //- Destructor
    virtual ~grass() = default;


    // Member Functions

    // Evaluation

        //- Return heat-transfer rate for speficied cells [J/s]
        virtual tmp<scalarField> Q(const labelList& cells) const;


    // I-O

        //- Read the dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace evapotranspirationHeatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
