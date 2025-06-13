/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024-2025 OpenCFD Ltd.
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
    Foam::functionObjects::bladeForces

Group
    grpForcesFunctionObjects

Description
    Computes forces and coefficients over a given list of patches by integrating
    pressure and viscous forces and moments.

    Forces and moments are output in their constituent components
    within the user-defined cylindrical coordinate system:
    - thrust
    - drag
    - torque

    Operands:
    \table
      Operand       | Type           | Location
      input         | -              | -
      output file   | dat            | postProcessing/\<FO\>/\<time\>/\<file\>s
    \endtable

    where \c \<file\>s:
    \verbatim
      force.dat        | Forces (thrust drag torque)
      C_d.dat          | Cd values (per bin)
      C_l.dat          | Cl values (per bin)
      C_p.dat          | Cp values (per bin)
    \endverbatim

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    <namePrefix>
    {
        // Mandatory entries
        type            bladeForces;
        libs            (forces);
        patches         (<wordRes>);

        // Optional entries
        writeFields     <bool>;
        fieldsInterval  <int>;
        useNamePrefix   <bool>;

        // Mandatory entries
        origin          (0 0 0);
        axis            (1 0 0);
        n               25;

        // Conditional optional entries
        p               <word>;
        U               <word>;
        rho             <word>;
        rhoInf          <scalar>;  // enabled if rho=rhoInf
        pRef            <scalar>;
        Uref            <scalar>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description               | Type | Reqd    | Deflt
      type       | Type name: bladeForces    | word | yes     | -
      libs       | Library name: forces      | word | yes     | -
      patches    | Names of operand patches  | wordRes | yes  | -
      writeFields | Flag to write surface with fields | bool | no | false
      fieldsInterval | Frequency of writeInterval for writeFields | int | no | 0
      useNamePrefix | Flag to include prefix for field names | bool | no | false
      outputName | Name for registered surface and VTP output | word | no | <name>
      origin  | Origin of cylindrical coordinate system | vector | yes | -
      axis    | Axis of cylindrical coordinate system   | vector | yes | -
      n       | Rotation speed [rev/sec]                | scalar | yes | -
      rpm     | Rotation speed [rev/min]                | scalar | no  | -
      p       | Name of pressure field      | word   | no   | p
      U       | Name of velocity field      | word   | no   | U
      rho     | Name of density field       | word   | no   | rho
      rhoInf  | Value of reference density  | scalar | cndtnl   | 1
      pRef    | Value of reference pressure | scalar | cndtnl   | 0
      Uref    | Magnitude of inlet reference axial velocity | scalar  | yes  | -

      radius   | The blade outer radius        | scalar | no  | 1
      nRadial  | Divisions in radial direction | label  | no  | 10
      lefthand | Using a left-hand blade         | bool | false  | false
      nearCellValue | Patch velocity extrapolated from fluid | bool | no | false
    \endtable

    Experimental entries (may be removed in the future):
    \table
      Property   | Description               | Type | Reqd    | Deflt
      geometricVelocity | Patch velocity based on position | bool | no | false
      mag.thrust | Ignore sign for thrust values | bool | false  | false
      mag.drag | Ignore sign for drag values   | bool | false  | false
    \endtable

    The inherited entries are elaborated in:
      - \link functionObject.H \endlink
      - \link writeFile.H \endlink
      - \link coordinateSystem.H \endlink

Note
  - For incompressible cases, set \c rho to \c rhoInf.
    You will then be required to provide a \c rhoInf
    value corresponding to the constant freestream density.
  - \c writeControl and \c writeInterval entries of function
    object do control when to output force/coefficients files and fields.

SourceFiles
    bladeForces.cxx

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_bladeForces_H
#define Foam_functionObjects_bladeForces_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "cylindricalCS.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class bladeForces Declaration
\*---------------------------------------------------------------------------*/

class bladeForces
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private Data

    // Results

        //- Sum of thrust forces (axial)
        scalar sumThrust_;

        //- Sum of drag forces (tangential)
        scalar sumDrag_;

        //- Sum of torque moments (tangential)
        scalar sumTorque_;

        //- Overall blade area (single side)
        scalar totalArea_;

        //- Overall drag coefficient (area-averaged of bin coefficients)
        scalar totalCd_;

        //- Overall lift coefficient (area-averaged of bin coefficients)
        scalar totalCl_;

        //- Overall pressure coefficient (area-averaged of bin coefficients)
        scalar totalCp_;

        //- Blade area (single side) per radial band
        scalarList bandArea_;

        //- Drag coefficient per radial band
        scalarList bandCd_;

        //- Lift coefficient per radial band
        scalarList bandCl_;

        //- Pressure coefficient per radial band
        scalarList bandCp_;


    // Geometric Information

        //- Selected operand patches
        labelList patchIDs_;

        //- Max (or reference) blade radius
        scalar refRadius_;

        //- Number of divisions in radial direction
        label nRadialDiv_;

        //- Rotational speed (revolutions per second)
        scalar revPerSec_;

        //- Cylindrical coordinate system used for force and torque
        coordSystem::cylindrical cylCoord_;


     // Read from dictionary

        //- Reference axial velocity
        scalar Uref_;

        //- Reference pressure
        scalar pRef_;

        //- Reference density (for incompressible)
        scalar rhoRef_;

        //- Name of pressure field (default: "p")
        word pName_;

        //- Name of velocity field (default: "U")
        word UName_;

        //- Name of density field (default: "rho")
        word rhoName_;

        //- Name for registred surface and VTP output.
        //  Default is the functionObject::name()
        word outputName_;

        //- Internal counter, incremented when write() is called.
        label writeCounter_;

        //- Write surface/fields interval.
        //  A value <= 1 means write surface/fields at each write().
        label fieldsInterval_;

        //- Flag of initialisation (internal)
        bool initialised_;

        //- Flag to write force and moment fields
        bool writeFields_;

        //- Using a left-hand blade (flips orientation of drag value)
        bool lefthand_;

        //- Extrapolate velocity to patch
        bool nearCellValue_;

        //- Blade speed based on position [experimental]
        bool useGeometricVelocity_;

        //- Ignore sign for thrust values [experimental]
        bool useMagThrust_;

        //- Ignore sign for drag values [experimental]
        bool useMagDrag_;


    // File streams

        //- File stream for forces, torque and summary
        autoPtr<OFstream> forceFilePtr_;

        //- File stream for areas (in radial bins)
        autoPtr<OFstream> areaFilePtr_;

        //- File stream for coefficient of drag (in radial bins)
        autoPtr<OFstream> CdFilePtr_;

        //- File stream for coefficient of lift (in radial bins)
        autoPtr<OFstream> ClFilePtr_;

        //- File stream for coefficient of pressure (in radial bins)
        autoPtr<OFstream> CpFilePtr_;


    // Private Member Functions

        //- Initialise containers and fields
        void initialise();

        //- Write surface (VTK format) and fields
        void writeSurface() const;


    // Evaluation

        //- Return the effective stress (viscous + turbulent) for patch
        tmp<symmTensorField> devRhoReff
        (
            const tensorField& gradUp,
            const label patchi
        ) const;

        //- Return dynamic viscosity field
        tmp<volScalarField> mu() const;

        //- Return rho if specified otherwise rhoRef
        tmp<volScalarField> rho() const;

        //- Return rho if specified otherwise rhoRef for patch
        tmp<scalarField> rho(const label patchi) const;

        //- Return rhoRef if the pressure field is dynamic (i.e. p/rho),
        //- otherwise return 1
        scalar rho(const volScalarField& p) const;


    // I-O

        //- Create the integrated-data files
        void createIntegratedDataFiles();

        //- Write integrated data to files
        void writeIntegratedDataFiles();


public:

    //- Runtime type information
    TypeName("bladeForces");


    // Constructors

        //- Construct from Time and dictionary
        bladeForces
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool readFields = true
        );

        //- Construct from objectRegistry and dictionary
        bladeForces
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool readFields = true
        );

        //- No copy construct
        bladeForces(const bladeForces&) = delete;

        //- No copy assignment
        void operator=(const bladeForces&) = delete;


    //- Destructor
    virtual ~bladeForces() = default;


    // Member Functions

        //- The integrated thrust force (axial)
        scalar thrust() const noexcept { return sumThrust_; }

        //- The integrated drag force (tangential)
        scalar drag() const noexcept { return sumDrag_; }

        //- The integrated torque
        scalar torque() const noexcept { return sumTorque_; }

        //- Calculate forces, torque, coefficients
        void calculate();


        //- Read the dictionary
        virtual bool read(const dictionary& dict);

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
