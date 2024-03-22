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
    Foam::evapotranspirationHeatTransferModel

Description
    Base class for evapotranspiration models
    to handle general evapotranspiration characteristics.

SourceFiles
    evapotranspirationHeatTransferModel.C
    evapotranspirationHeatTransferModelNew.C
    evapotranspirationHeatTransferModelTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_evapotranspirationHeatTransferModel_H
#define Foam_evapotranspirationHeatTransferModel_H

#include "dictionary.H"
#include "HashSet.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"
#include "OFstream.H"
#include "coordinateSystem.H"
#include "writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
            Class evapotranspirationHeatTransferModel Declaration
\*---------------------------------------------------------------------------*/

class evapotranspirationHeatTransferModel
{
    // Private Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Identifier of specified face zone
        label faceZoneId_;


protected:

    // Protected Member Functions

        //- Return requested field from the object registry
        //- or read+register the field to the object registry
        volScalarField& getOrReadField(const word& fieldName) const;

        //- Area-averaged heat flux due to short-wave
        //- solar rays on face zone [W/m^2]
        //  Short-wave: direct and diffusive solar rays
        scalar q() const;


public:

    //- Runtime type information
    TypeName("evapotranspirationHeatTransferModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            evapotranspirationHeatTransferModel,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (dict, mesh)
        );


    // Selectors

        //- Return a reference to the selected evapotranspiration model
        static autoPtr<evapotranspirationHeatTransferModel> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Constructors

        //- Construct from components
        evapotranspirationHeatTransferModel
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        evapotranspirationHeatTransferModel
        (
            const evapotranspirationHeatTransferModel&
        ) = delete;

        //- No copy assignment
        void operator=(const evapotranspirationHeatTransferModel&) = delete;


    //- Destructor
    virtual ~evapotranspirationHeatTransferModel() = default;


    // Member Functions

    // Access

        //- Return const reference to the mesh
        const fvMesh& mesh() const noexcept
        {
            return mesh_;
        }

        //- Return the face-zone identifier
        label faceZoneId() const noexcept
        {
            return faceZoneId_;
        }


    // Evaluation

        //- Return heat-transfer rate for speficied cells [J/s]
        virtual tmp<scalarField> Q(const labelList& cells) const = 0;


    // I-O

        //- Read the dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
