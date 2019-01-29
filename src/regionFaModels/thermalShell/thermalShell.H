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
    Foam::regionFaModels::thermalShell

Description
    Thermal shell finite area. It solves the energy Eq in 2D.
    The coupling with the 3D region is done through the temperatureCoupledBase,
    plus faOption is available to add extra sources on the shell such as
    external heat, joule heating, etc.


SourceFiles
    thermalShell.C
    thermalShellI.H

\*---------------------------------------------------------------------------*/

#ifndef thermalShell_H
#define thermalShell_H

#include "volFieldsFwd.H"
#include "thermalShellModel.H"
#include "solidProperties.H"
#include "faMesh.H"
#include "faCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{

/*---------------------------------------------------------------------------*\
                        Class thermalShell Declaration
\*---------------------------------------------------------------------------*/

class thermalShell
:
    public thermalShellModel
{
private:

    // Private member functions

        //- No copy construct
        thermalShell(const thermalShell&) = delete;

        //- No copy assignment
        void operator=(const thermalShell&) = delete;

        //- Initialize thermalShell
        void init();


protected:

    // Protected data

        // Solution parameters

            //- Number of non orthogonal correctors
            label nNonOrthCorr_;


        // Thermo properties

            //- Solid properties
            solidProperties thermo_;


        // Source term fields

            //- External surface energy source  / [J/m2/s]
            areaScalarField qs_;

            //- Temperature - wall [K]
            areaScalarField Tw1_;

            //- Thickness
            areaScalarField h_;


        // Sub models

            //- Pointer to radiation model
            //autoPtr<radiation::radiationModel> radiation_;


    // Protected member functions

        //- Heat transfe to primary region
        tmp<faScalarMatrix> q(areaScalarField& T);

        //- Read control parameters from dictionary
        virtual bool read(const dictionary& dict);


        // Equations

            //- Solve energy equation
            void solveEnergy();


public:

    //- Runtime type information
    TypeName("thermalShell");


    // Constructors


        //- Construct from components and dict
        thermalShell
        (
            const word& modelType,
            const fvMesh& mesh,
            const fvPatch& patch,
            const dictionary& dict
        );


    //- Destructor
    virtual ~thermalShell();


    // Member Functions


        // Fields

            //- Return the film specific heat capacity [J/kg/K]
            const tmp<areaScalarField> Cp() const;

            //- Return density [Kg/m3]
            const tmp<areaScalarField> rho() const;

            //- Return thermal conductivity [W/m/K]
            const tmp<areaScalarField> kappa() const;


        // Evolution

            //- Pre-evolve  thermal baffle
            virtual void preEvolveRegion();

            //- Evolve the thermal baffle
            virtual void evolveRegion();


       // I-O

            //- Provide some feedback
            virtual void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
