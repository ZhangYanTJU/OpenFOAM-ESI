/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenCFD Ltd.
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
    Foam::wallHeatFluxModels::gauge

Description
    This model computes the wall-heat flux at a selected patch using
    imaginary heat-flux gauges.

    Heat-flux gauges are measurement tools that are placed on walls. The
    temperature of the gauges is usually different from that of the surrounding
    wall - an order of magnitude lower. The model allows the user to
    specify the gauge temperature and other gauge properties such as
    absorptivity or emissivity without creating separate patches.
    The model calculates the net convective and radiative heat fluxes at the
    gauge locations, combines them, and writes the results to the output files.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    wallHeatFlux1
    {
        // Mandatory entries
        type        wallHeatFlux;
        libs        (fieldFunctionObjects);

        model       gauge;

        Tgauge      <scalar>;
        patch       <word>;

        // Optional entries
        absorptivity  <scalar>;
        emissivity    <scalar>;
        T             <word>;
        qin           <word>;
        alphat        <word>;
        convective    <bool>;
        radiative     <bool>;
        writeFields   <bool>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                           | Type | Reqd | Deflt
      type         | Type name: wallHeatFlux               | word |  yes  | -
      libs         | Library name: fieldFunctionObjects    | word |  yes  | -
      model      | Model name: gauge                     | word |  no  | wall
      Tgauge     | Gauge temperature (K)                  | scalar | yes | -
      patch      | Name of the patch to probe              | word | yes  | -
      absorptivity | Absorptivity of the gauge              | scalar | no | 1
      emissivity   | Emissivity of the gauge                | scalar | no | 1
      T           | Name of the temperature field          | word | no  | T
      qin | Name of the incident radiative heat flux field | word | no | qin
      alphat | Name of turbulent thermal diffusivity field | word | no | alphat
      convective  | Calculate convective heat flux         | bool | no | false
      radiative   | Calculate radiative heat flux          | bool | no | false
      writeFields | Write the fields to file                | bool | no | false
    \endtable

    The inherited entries are elaborated in:
     - \link fvMeshFunctionObject.H \endlink
     - \link writeFile.H \endlink
     - \link pointProberBase.H \endlink
     - \link patchProber.H \endlink

SourceFiles
    gauge.C

\*---------------------------------------------------------------------------*/

#ifndef wallHeatFluxModels_gauge_H
#define wallHeatFluxModels_gauge_H

#include "wallHeatFluxModel.H"
#include "writeFile.H"
#include "patchProber.H"
#include "wallDist.H"
#include "turbulentFluidThermoModel.H"
#include "scalarIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallHeatFluxModels
{

/*---------------------------------------------------------------------------*\
                            Class gauge Declaration
\*---------------------------------------------------------------------------*/

class gauge
:
    public wallHeatFluxModel,
    public patchProber
{
    // Private Data

        //- Const reference to the thermo database
        const fluidThermo& thermo_;

        //- Net convective gauge heat flux at patch-probe locations
        std::unique_ptr<scalarIOField> qConvPtr_;

        //- Net radiative gauge heat flux at patch-probe locations
        std::unique_ptr<scalarIOField> qRadPtr_;

        //- Flag to output the net convective heat flux data
        autoPtr<OFstream> convectiveFilePtr_;

        //- Flag to output the net radiative heat flux data
        autoPtr<OFstream> radiativeFilePtr_;

        //- Indices of the patch cells that will be probed
        labelField cells_;

        //- Identifier of the faces that are actively sampled
        boolList activeFaces_;

        //- Name of the incident radiative heat-flux field
        word qinName_;

        //- Name of the temperature field
        word Tname_;

        //- Name of the turbulent thermal diffusivity field
        word alphatName_;

        //- Common and constant temperature of the gauges
        scalar Tgauge_;

        //- Radiant exitance
        scalar Me_;

        //- Absorptivity of the gauge surfaces
        scalar a_;

        //- Emissivity of the gauge surfaces
        scalar e_;

        //- Index of the operand patch
        label patchID_;

        //- Number of specified patch probes
        label szProbes_;

        //- Flag to output the net convective/radiative heat flux fields
        bool writeFields_;


    // Private Member Functions

        //- Return the radiant-exitance term of the radiative heat flux
        scalar calcRadiantExitance() const;

        //- Return the indices of the cells that contain the probed patch faces
        tmp<labelField> identifyProbeCells() const;

        //- Write the file-header information
        bool writeFileHeader(Ostream& os);

        //- Calculate the net convective heat flux
        bool calcConvectiveHeatFlux();

        //- Return the dT/dn term of the convective heat flux
        tmp<scalarField> calcdTdn() const;

        //- Calculate the net radiative heat flux
        bool calcRadiativeHeatFlux();


public:

    //- Runtime type information
    TypeName("gauge");


    // Constructors

        //- Construct from components
        gauge
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& name,
            const word objName,
            functionObjects::stateFunctionObject& state
        );


    //- Destructor
    virtual ~gauge() = default;


    // Member Functions

    // Evaluation

        //- Read the settings
        virtual bool read(const dictionary& dict);

        //- Calculate the heat-flux data
        virtual bool execute();

        //- Write the heat-flux data
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wallHeatFluxModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
