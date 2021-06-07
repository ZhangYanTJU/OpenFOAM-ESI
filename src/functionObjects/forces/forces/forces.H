/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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
    Foam::functionObjects::forces

Group
    grpForcesFunctionObjects

Description
    Calculates the forces and moments by integrating the pressure and
    skin-friction forces over a given list of patches, and the resistance
    from porous zones.

    Forces and moments are calculated in a global Cartesian coordinate system
    by default, or using a user-defined system.  Contributions can be 'binned'
    according to a user-defined number of uniform-width collection zones (bins)
    that span the input geometries, oriented by a user-defined direction vector.

    Results are written to multiple files as a function of time in the
    postProcessing/\<functionObjectName\> directory:
    - force.dat          : forces
    - moment.dat         : moments
    - forceBin.dat       : force bins
    - momentBin.dat      : moment bins

Usage
    Example of function object specification:
    \verbatim
    forces1
    {
        type        forces;
        libs        (forces);
        ...
        log         yes;
        writeFields yes;
        patches     (walls);

        binData
        {
            nBin        20;
            direction   (1 0 0);
            cumulative  yes;
        }
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: forces       | yes         |
        log          | Write force data to standard output | no | no
        writeFields  | Write the force and moment fields | no | no
        patches      | Patches included in the forces calculation | yes |
        p            | Pressure field name     | no          | p
        U            | Velocity field name     | no          | U
        rho          | Density field name (see below) | no   | rho
        CofR         | Centre of rotation (see below) | no   |
        porosity     | flag to include porosity contributions | no | no
        directForceDensity | Force density supplied directly (see below)|no|no
        fD           | Name of force density field (see below) | no | fD
    \endtable

    Bin data is optional, but if the dictionary is present, the entries must
    be defined according o
    \table
        nBin         | number of data bins     | yes         |
        direction    | direction along which bins are defined | yes |
        cumulative   | bin data accumulated with incresing distance | yes |
    \endtable

Note
  - For incompressible cases, set \c rho to \c rhoInf.  You will then be
    required to provide a \c rhoInf value corresponding to the free-stream
    constant density.
  - If the force density is supplied directly, set the \c directForceDensity
    flag to 'yes', and supply the force density field using the \c
    fDName entry
  - The centre of rotation (CofR) for moment calculations can either be
    specified by an \c CofR entry, or be taken from origin of the local
    coordinate system.  For example,
    \verbatim
        CofR        (0 0 0);
    \endverbatim
    or
    \verbatim
        origin  (0 0 0);
        e1      (0 1 0);
        e3      (0 0 1);
    \endverbatim
    or
    \verbatim
        coordinateSystem
        {
            origin  (0 0 0);
            rotation
            {
                type    axes;
                e3      (0 0 1);
                e1      (1 0 0);
            }
        }
    \endverbatim

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    forces.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_forces_H
#define functionObjects_forces_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "cartesianCS.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class forces Declaration
\*---------------------------------------------------------------------------*/

class forces
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected data

        //- Pressure, viscous and porous force per bin
        List<Field<vector>> forces_;

        //- Pressure, viscous and porous moment per bin
        List<Field<vector>> moments_;

        //- Integrated force constituents
        //  [0] = total; [1] = pressure; [2] = viscous; (optional: [3] = porous)
        List<vector> sumForces_;

        //- Integrated moment constituents
        //  [0] = total; [1] = pressure; [2] = viscous; (optional: [3] = porous)
        List<vector> sumMoments_;


        // File streams

            //- Forces
            autoPtr<OFstream> forceFilePtr_;

            //- Moments
            autoPtr<OFstream> momentFilePtr_;


        // Read from dictionary

            //- Coordinate system used when evaluating forces/moments
            coordSystem::cartesian coordSys_;

            //- Patches to integrate forces over
            labelHashSet patchSet_;

            //- Reference density needed for incompressible calculations
            scalar rhoRef_;

            //- Reference pressure
            scalar pRef_;

            //- Name of pressure field
            word pName_;

            //- Name of velocity field
            word UName_;

            //- Name of density field (optional)
            word rhoName_;

            //- The name of the force density (fD) field
            word fDName_;

            //- Is the force density being supplied directly?
            bool directForceDensity_;

            //- Flag to include porosity effects
            bool porosity_;

            //- Write fields flag
            bool writeFields_;

            //- Initialised flag
            bool initialised_;


    // Protected Member Functions

        //- Create a field name
        word fieldName(const word& name) const;

        //- Set the co-ordinate system from dictionary and axes names
        void setCoordinateSystem
        (
            const dictionary& dict,
            const word& e3Name = word::null,
            const word& e1Name = word::null
        );

        //- Initialise the fields
        void initialise();

        //- Reset containers and fields
        void reset();


        // Evaluation

            //- Return the effective stress (viscous + turbulent)
            tmp<volSymmTensorField> devRhoReff() const;

            //- Return dynamic viscosity field
            tmp<volScalarField> mu() const;

            //- Return rho if specified otherwise rhoRef
            tmp<volScalarField> rho() const;

            //- Return rhoRef if the pressure field is
            //- dynamic (i.e. p/rho), otherwise return 1
            scalar rho(const volScalarField& p) const;

            void applyBins
            (
                const vectorField& Md,
                const vectorField& fN,
                const vectorField& fT,
                const vectorField& fP,
                const vectorField& d
            );

            //- Add patch contributions to force and moment fields
            void addToFields
            (
                const label patchi,
                const vectorField& Md,
                const vectorField& fN,
                const vectorField& fT
            );

            //- Add cell contributions to force and
            //- moment fields, and include porosity effects
            void addToFields
            (
                const labelList& cellIDs,
                const vectorField& Md,
                const vectorField& fN,
                const vectorField& fT,
                const vectorField& fP
            );

            //- Return integrated {total, pressure, viscous, porous} components
            List<vector> integrateData
            (
                const vectorField& fm0,
                const vectorField& fm1,
                const vectorField& fm2
            ) const;


        // Write

            //- Create the integrated-data files
            bool createDataFiles();

            //- Write header for an integrated-data file
            bool writeDataFileHeader
            (
                const word& header,
                OFstream& os
            ) const;

            //- Write integrated data to files
            bool writeDataFiles();

            //- Write integrated data to a file
            bool writeDataFile
            (
                const List<vector>& fm,
                autoPtr<OFstream>& osPtr
            ) const;

            //- Write integrated data to stream
            void logData
            (
                const string& descriptor,
                const List<vector>& fm
            ) const;


    // Generated Methods

        //- No copy construct
        forces(const forces&) = delete;

        //- No copy assignment
        void operator=(const forces&) = delete;


public:

    //- Runtime type information
    TypeName("forces");


    // Constructors

        //- Construct from Time and dictionary
        forces
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool readFields = true
        );

        //- Construct from objectRegistry and dictionary
        forces
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool readFields = true
        );


    //- Destructor
    virtual ~forces() = default;


    // Member Functions

        //- Read the dictionary
        virtual bool read(const dictionary& dict);

        //- Calculate forces and moments
        virtual void calcForcesMoments();

        //- Return the total force
        virtual vector forceEff() const;

        //- Return the total moment
        virtual vector momentEff() const;

        //- Execute the function object
        virtual bool execute();

        //- Write to data files/fields and to streams
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
