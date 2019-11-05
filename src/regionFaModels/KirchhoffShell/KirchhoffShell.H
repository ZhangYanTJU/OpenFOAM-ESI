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
    Foam::regionFaModels::KirchhoffShell

Description
    Foam::regionFaModels::KirchhoffShell


SourceFiles
    KirchhoffShell.C

\*---------------------------------------------------------------------------*/

#ifndef KirchhoffShell_H
#define KirchhoffShell_H

#include "volFieldsFwd.H"
#include "vibrationShellModel.H"
#include "faMesh.H"
#include "faCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{

/*---------------------------------------------------------------------------*\
                        Class KirchhoffShell Declaration
\*---------------------------------------------------------------------------*/

class KirchhoffShell
:
    public vibrationShellModel
{
private:

    // Private data
    
        //- Damping coefficients [1/s]
        dimensionedScalar f0_;
        dimensionedScalar f1_;
        dimensionedScalar f2_;
        
        
    // Private member functions

        //- No copy construct
        KirchhoffShell(const KirchhoffShell&) = delete;

        //- No copy assignment
        void operator=(const KirchhoffShell&) = delete;

        //- Initialize KirchhoffShell
        void init();


protected:

    // Protected data

        // Solution parameters

            //- Number of non orthogonal correctors
            label nNonOrthCorr_;


        // Source term fields

            //- External surface source [Pa]
            areaScalarField ps_;

            //- Thickness [m]
            areaScalarField h_;
            
            //- Cached Laplacian of w for dumping model
            areaScalarField laplaceW_;
            
            //- Cached Laplacian of  laplacian w for dumping model
            areaScalarField laplace2W_;



    // Protected member functions


        //- Read control parameters from dictionary
        virtual bool read(const dictionary& dict);


        // Equations

            //- Solve energy equation
            void solveDisplacement();


public:

    //- Runtime type information
    TypeName("KirchhoffShell");


    // Constructors


        //- Construct from components and dict
        KirchhoffShell
        (
            const word& modelType,
            const fvPatch& patch,
            const dictionary& dict
        );


    //- Destructor
    virtual ~KirchhoffShell();


    // Member Functions


        // Fields

            //- Return stiffness
            const tmp<areaScalarField> D() const;

            //- Return density [Kg/m3]
            const tmp<areaScalarField> rho() const;


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
