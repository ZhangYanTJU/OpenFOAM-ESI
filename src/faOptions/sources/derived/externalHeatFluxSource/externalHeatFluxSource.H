/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
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
    Foam::fa::externalHeatFluxSource


Description
    This faOption applies a heat flux condition for a faMesh region
    to temperature on an external wall in one of three modes:

      - fixed power: supply Q
      - fixed heat flux: supply q
      - fixed heat transfer coefficient: supply h and Ta

    where:
    \vartable
        Q  | Power [W]
        q  | Heat flux [W/m^2]
        h  | Heat transfer coefficient [W/m^2/K]
        Ta | Ambient temperature [K]
    \endvartable

    The ambient temperature Ta is specified as a Foam::Function1 of time but
    uniform in space.

Usage
    \table
    Property    | Description                              | Required | Default
    mode        | 'power', 'flux' or 'coefficient'         | yes |
    Q           | Power [W]                                | for mode 'power' |
    q           | Heat flux [W/m^2]                        | for mode 'flux' |
    h           | Heat transfer coefficient [W/m^2/K] | for mode 'coefficient' |
    Ta          | Ambient temperature [K]             | for mode 'coefficient' |
    relaxation  | Relaxation for the wall temperature      | no  | 1
    emissivity  | Surface emissivity for radiative flux to ambient | no | 0
    region      | 2D region to be applied                  | yes
    T           | Temperature field name                   | yes
    \endtable

    Example of the source specification:
    \verbatim
    <patchName>
    {
        type                                externalHeatFluxSource;
        active                              true;
        selectionMode                       all;

        region                              ceilingShell;
        T                                   Ts_ceilingShell;

        mode                                coefficient;
        Ta                                  constant 280;
        h                                   10;
    }
    \endverbatim


SourceFiles
    externalHeatFluxSource.C

\*---------------------------------------------------------------------------*/

#ifndef fa_externalHeatFluxSource_H
#define fa_externalHeatFluxSource_H

#include "faOption.H"
#include "Function1.H"
#include "areaFields.H"
#include "faceSetOption.H"
#include "faCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fa
{

/*---------------------------------------------------------------------------*\
                    Class externalHeatFluxSource Declaration
\*---------------------------------------------------------------------------*/

class externalHeatFluxSource
:
    public faceSetOption
{
public:

    // Public data

        //- Operation mode enumeration
        enum operationMode
        {
            fixedPower,                 //!< Fixed heat power [W]
            fixedHeatFlux,              //!< Fixed heat flux [W/m2]
            fixedHeatTransferCoeff,     //!< Fixed heat transfer coefficient
        };

        static const Enum<operationMode> operationModeNames;


private:

    // Private data

        //- Operation mode
        operationMode mode_;

        //- Name of temperature field - default = "T" (optional)
        word TName_;

        //- Heat power [W]
        scalar Q_;

        //- External Heat flux [W/m2]
        scalar q_;

        //- Heat transfer coefficient [W/m2K]
        scalar h_;

        //- Ambient temperature [K]
        autoPtr<Function1<scalar>> Ta_;

        //- Optional surface emissivity for radiative transfer to ambient
        scalar emissivity_;


public:

    //- Runtime type information
    TypeName("externalHeatFluxSource");


    // Constructors

         //- Construct from explicit source name and mesh
        externalHeatFluxSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh,
            const fvPatch& patch
        );



    //- Destructor
    virtual ~externalHeatFluxSource();


    // Member Functions

        // Evaluate

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const areaScalarField& h,
                const areaScalarField& rho,
                faMatrix<scalar>& eqn,
                const label fieldi
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
