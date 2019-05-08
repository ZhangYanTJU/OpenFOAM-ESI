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
    Foam::regionModels::thermalShellModels::thermalShellModel

Description
    Base class for thermal 2D shells

SourceFiles
    thermalShellModel.C

\*---------------------------------------------------------------------------*/

#ifndef thermalShellModel_H
#define thermalShellModel_H

#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "areaFieldsFwd.H"
#include "volFieldsFwd.H"
#include "regionFaModel.H"
#include "faOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//using Foam::fa;

namespace Foam
{
namespace regionModels
{

/*---------------------------------------------------------------------------*\
                      Class thermalShellModel Declaration
\*---------------------------------------------------------------------------*/

class thermalShellModel
:
    public regionFaModel
{
private:

    // Private Member Functions

        //- No copy construct
        thermalShellModel(const thermalShellModel&) = delete;

        //- No copy assignment
        void operator=(const thermalShellModel&) = delete;

        //- Initialize thermal Baffle
        void init();


protected:

    // Protected Data

        //- Name of the temperature field
        word TName_;

        //- Primary region temperature
        const volScalarField& Tp_;

        //- Shell temperature
        areaScalarField T_;

        //- Pointer to faOptions
        Foam::fa::options& faOptions_;


    // Protected Member Functions

        //- Read control parameters from dictionary
        virtual bool read(const dictionary&);


public:

    //- Runtime type information
    TypeName("thermalShellModel");


    // Declare runtime constructor selection tables

         declareRunTimeSelectionTable
         (
             autoPtr,
             thermalShellModel,
             dictionary,
             (
                const word& modelType,
                const fvMesh& mesh,
                const fvPatch& patch,
                const dictionary& dict
             ),
             (modelType, mesh, patch, dict)
         );


    // Constructors


        //- Construct from type name and mesh and dict
        thermalShellModel
        (
            const word& modelType,
            const fvMesh& mesh,
            const fvPatch& patch,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected model using dictionary
        static autoPtr<thermalShellModel> New
        (
            const fvMesh& mesh,
            const fvPatch& patch,
            const dictionary& dict
        );


    //- Destructor
    virtual ~thermalShellModel();


    // Member Functions


        // Access functions

            //- Return primary region temperature
            const volScalarField& Tp() const;

            //- Return shell temperature
            const areaScalarField& T() const;

            //- Return faOptions
            Foam::fa::options& faOptions();


        // Evolution

            //- Pre-evolve region
            virtual void preEvolveRegion();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
