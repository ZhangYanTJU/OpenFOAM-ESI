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
    Foam::fv::evapotranspirationHeatTransfer

Description
    Applies sources on temperature (\c T - incompressible)
    or energy (\c h/e - compressible) equation to incorporate
    evapotranspiration heat-transfer effects from the specified
    plant canopy. Heat transfer is usually calculated based on
    empirical relations between plants and solar radiation.

Usage
    Example by using \c constant/fvOptions:
    \verbatim
    evapotranspirationHeatTransfer1
    {
        // Mandatory entries
        type                evapotranspirationHeatTransfer;
        model               <word>;

        // Model-specific entries
        ...

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                         | Type | Reqd | Deflt
      type     | Type name: evapotranspirationHeatTransfer  | word | yes  | -
      model    | Name of heat-transfer model        | word  | yes  | -
    \endtable

    Options for the \c model entry:
    \verbatim
      tree   | Regression-based heat-transfer model for trees
      grass  | Heat-transfer model for grass
    \endverbatim

    The inherited entries are elaborated in:
      - \link fvOption.H \endlink
      - \link cellSetOption.H \endlink

Note
  - Evapotranspiration mass transfer is not modelled.

SourceFiles
    evapotranspirationHeatTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fv_evapotranspirationHeatTransfer_H
#define Foam_fv_evapotranspirationHeatTransfer_H

#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class evapotranspirationHeatTransferModel;

namespace fv
{

/*---------------------------------------------------------------------------*\
                Class evapotranspirationHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class evapotranspirationHeatTransfer
:
    public fv::cellSetOption
{
    // Private Data

        //- Runtime-selectable evapotranspiration heat transfer model
        autoPtr<evapotranspirationHeatTransferModel> ethtModelPtr_;


public:

    //- Runtime type information
    TypeName("evapotranspirationHeatTransfer");


    // Constructors

        //- Construct from explicit source name and mesh
        evapotranspirationHeatTransfer
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        evapotranspirationHeatTransfer(const evapotranspirationHeatTransfer&) =
            delete;

        //- No copy assignment
        void operator=(const evapotranspirationHeatTransfer&) = delete;


    //- Destructor
    virtual ~evapotranspirationHeatTransfer();


    // Member Functions

        //- Add explicit contribution to energy equation
        //- for incompressible flow computations
        virtual void addSup
        (
            fvScalarMatrix& eqn,
            const label fieldi
        );

        //- Add explicit contribution to energy equation
        //- for compressible flow computations
        virtual void addSup
        (
            const volScalarField& rho,
            fvScalarMatrix& eqn,
            const label fieldi
        );

        //- Add explicit contribution to energy equation
        //- for multiphase flow computations
        virtual void addSup
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            fvScalarMatrix& eqn,
            const label fieldi
        )
        {
            NotImplemented;
        }

        //- Read dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
