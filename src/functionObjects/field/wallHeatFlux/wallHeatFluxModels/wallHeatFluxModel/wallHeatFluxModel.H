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

Namespace
    Foam::wallHeatFluxModels

Description
    A namespace for various heat-flux model implementations.

Class
    Foam::wallHeatFluxModel

Description
    A base class for heat-flux models.

SourceFiles
    wallHeatFluxModel.C
    wallHeatFluxModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_wallHeatFluxModel_H
#define Foam_wallHeatFluxModel_H

#include "writeFile.H"
#include "volFieldsFwd.H"
#include "HashSet.H"
#include "stateFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                   Class wallHeatFluxModel Declaration
\*---------------------------------------------------------------------------*/

class wallHeatFluxModel
:
    public functionObjects::writeFile
{
    // Private Data

        //- Const reference to the mesh
        const fvMesh& mesh_;

        //- Reference to the state function object
        functionObjects::stateFunctionObject& state_;

        //- Name of the function object
        const word objName_;


public:

    //- Runtime type information
    TypeName("wallHeatFluxModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            wallHeatFluxModel,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh,
                const word& name,
                const word objName,
                functionObjects::stateFunctionObject& state
            ),
            (dict, mesh, name, objName, state)
        );


    // Selectors

        //- Return a reference to the selected heat-flux model
        static autoPtr<wallHeatFluxModel> New
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& name,
            const word objName,
            functionObjects::stateFunctionObject& state
        );


    // Generated Methods

        //- No copy construct
        wallHeatFluxModel(const wallHeatFluxModel&) = delete;

        //- No copy assignment
        void operator=(const wallHeatFluxModel&) = delete;


    // Constructors

        //- Construct from components
        wallHeatFluxModel
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& name,
            const word objName,
            functionObjects::stateFunctionObject& state
        );


    //- Destructor
    virtual ~wallHeatFluxModel() = default;


    // Member Functions

    // Access

        //- Return const reference to the mesh
        const fvMesh& mesh() const noexcept { return mesh_; }

        //- Return const reference to the state function object
        functionObjects::stateFunctionObject& state() const noexcept
        {
            return state_;
        }

        //- Return const reference to the function-object name
        const word& objName() const noexcept { return objName_; }


    // Evaluation

        //- Read the settings
        virtual bool read(const dictionary& dict) = 0;

        //- Calculate the heat-flux data
        virtual bool execute() = 0;

        //- Write the heat-flux data
        virtual bool write() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
