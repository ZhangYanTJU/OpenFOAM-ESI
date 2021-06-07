/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::functionObjects::forceCoeffs

Group
    grpForcesFunctionObjects

Description
    Computes force and moment coefficients over a given
    list of patches, and optionally over given porous zones.
    The following coefficients can be selected and output:
    \verbatim
      Cd      | Drag coefficient
      Cs      | Side-force coefficient
      Cl      | Lift coefficient
      CmRoll  | Roll-moment coefficient
      CmPitch | Pitch-moment coefficient
      CmYaw   | Yaw-moment coefficient
    \endverbatim

    The force coefficients can also be optionally output
    in terms of their front and rear axle constituents:
    \verbatim
        Cd{f,r} = 0.5*Cd {+,-} CmRoll
        Cs{f,r} = 0.5*Cs {+,-} CmYaw
        Cl{f,r} = 0.5*Cl {+,-} CmPitch
    \endverbatim
    where \c f and \c r represent front and rear axles, respectively.

    Force and moment coefficients are output
    in their total and constituent components:
    - total force and moment coefficients
    - pressure contributions
    - viscous contributions
    - porous resistance contributions (optional)

    Force and moment coefficients can be computed and output in:
    - the global Cartesian coordinate system (default)
    - a user-defined Cartesian coordinate system

    Force and moment coefficients can also be equidistant-spatially
    'binned' over given patches. See the following link for the details:
     - \link binModel.H \endlink

    Operands:
    \table
      Operand       | Type           | Location
      input         | -              | -
      output file   | dat | $FOAM_CASE/postProcessing/\<FO\>/\<time\>/\<file\>s
      output field  | volVectorField | $FOAM_CASE/\<time\>/\<outField\>s
    \endtable

    where \c \<file\>s:
    \verbatim
      coefficient.dat  | Integrated coefficients over all patches
    \endverbatim

    If a bin model is activated, additional \c \<file\>s are output:
    \verbatim
      CdBin.dat        | Drag coefficient bins
      CsBin.dat        | Side-force coefficient bins
      ClBin.dat        | Lift coefficient bins
      CmRollBin.dat    | Roll-moment coefficient bins
      CmPitchBin.dat   | Pitch moment coefficient bins
      CmYawBin.dat     | Yaw-moment coefficient bins
    \endverbatim

    where \c \<outField\>s:
    \verbatim
      forceCoeff       | Force coefficient field
      momentCoeff      | Moment coefficient field
    \endverbatim

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    forceCoeffs1
    {
        // Mandatory entries (unmodifiable)
        type                forceCoeffs;
        libs                (forces);

        // Mandatory entries (runtime modifiable)
        patches             (<patch1> ... <patchN>); // (wall1 "(wall2|wall3)");
        magUInf             100;
        lRef                3.5;
        Aref                2.2;

        // Optional entries (runtime modifiable)
        directForceDensity  false;
        porosity            false;
        writeFields         false;

        forceCoeffs         (Cd Cs Cl);
        momentCoeffs        (CmRoll CmPitch CmYaw);
        frontRearForceCoeffs (Cdf Csf Clf Cdr Csr Clr);

        // Conditional mandatory entries (runtime modifiable)

            // Cartesian coordinate system specification when evaluating
            // force and moment coefficients, either of the below

            // Option-1, i.e. the centre of rotation
            // by inherently using e3=(0 0 1) and e1=(1 0 0)
            CofR                (0 0 0); // Centre of rotation
            dragDir             (1 0 0);
            liftDir             (0 0 1);

            // Option-2, i.e. local coordinate system specification
            origin              (0 0 0);
            e1                  (1 0 0);
            e3                  (0 0 1); // (e1, e2) or (e2, e3) or (e3, e1)

            // Option-3, i.e. general coordinate system specification
            coordinateSystem
            {
                type            cartesian;
                origin          (0 0 0);
                rotation
                {
                    type        axes;
                    e3          (0 0 1);
                    e1          (1 0 0); // (e1, e2) or (e2, e3) or (e3, e1)
                }
            }

        // Conditional optional entries (runtime modifiable)

            // if directForceDensity == true
            fD              <fDName>;

            // if directForceDensity == false
            p               <pName>;
            U               <UName>;
            rho             <rhoName>;
            rhoInf          1; // redundant for incompressible-flow cases
            pRef            0;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description               | Type     | Reqd | Deflt
      type       | Type name: forceCoeffs    | word     | yes  | -
      libs       | Library name: forces      | word     | yes  | -
      patches    | Names of operand patches  | wordList | yes  | -
      magUInf    | Reference velocity magnitude       | scalar | yes | -
      lRef       | Reference length scale for moment          <!--
                 -->  calculations           | scalar | yes    | -
      Aref       | Reference area            | scalar | yes    | -
      directForceDensity | Flag to directly supply force      <!--
                 --> density - see below     | bool   | no     | false
      porosity   | Flag to include porosity contributions | bool | no | false
      writeFields | Flag to write force and moment fields | bool | no | false
      forceCoeffs  | Names of operand force coefficients      <!--
                 --> - see below             | wordList | no   | all
      momentCoeffs | Names of operand moment coefficients     <!--
                 --> - see below             | wordList | no   | all
      frontRearForceCoeffs | Names of operand front- and      <!--
                 --> rear-axle force coefficients - see below <!--
                 -->                         | wordList | no   | none
      CofR    | Centre of rotation           | vector | cndtnl | -
      origin  | Origin of coordinate system  | vector | cndtnl | -
      e3      | e3 coordinate axis           | vector | cndtnl | -
      e1      | e1 coordinate axis           | vector | cndtnl | -
      coordinateSystem | Coordinate system specifier | dictionary | cndtnl | -
      fD      | Name of force density field  | word   | cndtnl-no | fD
      p       | Name of pressure field       | word   | cndtnl-no | p
      U       | Name of velocity field       | word   | cndtnl-no | U
      rho     | Name of density field        | word   | cndtnl-no | rho
      rhoInf  | Value of reference density   | scalar | cndtnl-yes | -
      pRef    | Value of reference pressure  | scalar | cndtnl-no | 0
    \endtable

    Options for the \c forceCoeffs entry:
    \verbatim
      Cd      | Drag coefficient
      Cs      | Side-force coefficient
      Cl      | Lift coefficient
    \endverbatim

    Options for the \c momentCoeffs entry:
    \verbatim
      CmRoll  | Roll-moment coefficient
      CmPitch | Pitch-moment coefficient
      CmYaw   | Yaw-moment coefficient
    \endverbatim

    Options for the \c frontRearForceCoeffs entry:
    \verbatim
      Cdf     | Front-axle drag coefficient
      Csf     | Front-axle side-force coefficient
      Clf     | Front-axle lift coefficient
      Cdr     | Rear-axle drag coefficient
      Csr     | Rear-axle side-force coefficient
      Clr     | Rear-axle lift coefficient
    \endverbatim

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link writeFile.H \endlink
     - \link binModel.H \endlink
     - \link coordinateSystem.H \endlink

Note
  - \c rhoInf is always redundant for incompressible computations.
    That is, \c rhoInf is always equal to 1 in incompressible
    computations no matter which input value is assigned to \c rhoInf.
    The value of \c rhoInf is only used for compressible computations.
  - \c writeControl and \c writeInterval entries of function
    object do control when to output force and moment files.
  - \c writeControl and \c writeInterval entries of \c controlDict
    do control when to output force and moment fields.
  - Binning is not supported for force and moment fields,
    and for front-axle and rear-axle coefficients.
  - Input for force/moment coefficient directions
    require an origin and two orthogonal directions where
    the remaining orthogonal direction is automatically computed.
  - The default direction relations are shown below:

    \table
      Property     | Description           | Alias | Direction
      dragDir      | Drag direction        | e1    | (1 0 0)
      sideDir      | Side force direction  | e2    | (0 1 0)
      liftDir      | Lift direction        | e3    | (0 0 1)
      rollAxis     | Roll axis             | e1    | (1 0 0)
      pitchAxis    | Pitch axis            | e2    | (0 1 0)
      yawAxis      | Yaw axis              | e3    | (0 0 1)
    \endtable

See also
  - Foam::functionObject
  - Foam::functionObjects::forces
  - Foam::functionObjects::writeFile
  - Foam::binModel

SourceFiles
    forceCoeffs.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_forceCoeffs_H
#define functionObjects_forceCoeffs_H

#include "forces.H"
#include "HashSet.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class forceCoeffs Declaration
\*---------------------------------------------------------------------------*/

class forceCoeffs
:
    public forces
{
    // Private Enumerations

        //- Options for the force coefficients
        enum forceCoeffType
        {
            CD,        //!< Drag coefficient
            CS,        //!< Side-force coefficient
            CL         //!< Lift coefficient
        };

        //- Names for forceCoeffType
        static const Enum<forceCoeffType> forceCoeffTypeNames_;

        //- Options for the moment coefficients
        enum momentCoeffType
        {
            CMROLL,    //!< Roll-moment coefficient
            CMPITCH,   //!< Pitch-moment coefficient
            CMYAW,     //!< Yaw-moment coefficient
        };

        //- Names for momentCoeffType
        static const Enum<momentCoeffType> momentCoeffTypeNames_;

        //- Options for the front- and rear-axle force coefficients
        enum frontRearForceCoeffType
        {
            CDF,       //!< Front-axle drag coefficient
            CSF,       //!< Front-axle side-force coefficient
            CLF,       //!< Front-axle lift coefficient
            CDR,       //!< Rear-axle drag coefficient
            CSR,       //!< Rear-axle side-force coefficient
            CLR        //!< Rear-axle lift coefficient
        };

        //- Names for frontRearForceCoeffType
        static const Enum<frontRearForceCoeffType>
            frontRearForceCoeffTypeNames_;


    // Private Data

        //- Force coefficients per bin (i.e. Cd, Cs, Cl)
        //  [0] = pressure; [1] = viscous; (optional: [2] = porous)
        //  e.g. CFs_[0][1] is Cd.viscous bins
        List<List<scalarField>> CFs_;

        //- Moment coefficients per bin (i.e. CmRoll, CmPitch, CmYaw)
        //  [0] = pressure; [1] = viscous; [2] = porous
        //  e.g. CMs_[2][0] is CmYaw.pressure bins
        List<List<scalarField>> CMs_;

        //- Integrated force coefficients (i.e. Cd, Cs, Cl)
        //  [0] = total; [1] = pressure; [2] = viscous; (optional: [3] = porous)
        //  e.g. sumCFs_[0][2] is Cd.viscous
        List<List<scalar>> sumCFs_;

        //- Integrated moment coefficients (i.e. CmRoll, CmPitch, CmYaw)
        //  [0] = total; [1] = pressure; [2] = viscous; (optional: [3] = porous)
        //  e.g. sumCMs_[0][1] is CmRoll.pressure
        List<List<scalar>> sumCMs_;

        //- Integrated front-axle force coefficients (i.e. Cdf Csf Clf)
        //  [0] = total; [1] = pressure; [2] = viscous; (optional: [3] = porous)
        //  e.g. sumFrontCFs_[2][0] is Clf.total
        List<List<scalar>> sumFrontCFs_;

        //- Integrated rear-axle force coefficients (i.e. Cdr Csr Clr)
        //  [0] = total; [1] = pressure; [2] = viscous; (optional: [3] = porous)
        //  e.g. sumRearCFs_[2][2] is Clr.viscous
        List<List<scalar>> sumRearCFs_;

        //- Output force coefficients
        wordHashSet forceCoeffNames_;

        //- Output moment coefficients
        wordHashSet momentCoeffNames_;

        //- Output front and rear force coefficients
        wordHashSet frontRearForceCoeffNames_;

        //- Identifiers for operand force coefficients
        std::vector<label> idCFs_;

        //- Identifiers for operand moment coefficients
        std::vector<label> idCMs_;


        // File streams

            //- File stream for integrated operand coefficients
            autoPtr<OFstream> coeffFilePtr_;


        // Reference scales

            //- Reference velocity magnitude [m/s]
            scalar magUInf_;

            //- Reference length scale [m]
            scalar lRef_;

            //- Reference area [m^2]
            scalar Aref_;


        //- Flag to log front-axle and rear-axle coefficients (internal)
        bool calcFrontRear_;

        //- Flag of initialisation (internal)
        bool initialised_;


    // Generated Methods

        //- No copy construct
        forceCoeffs(const forceCoeffs&) = delete;

        //- No copy assignment
        void operator=(const forceCoeffs&) = delete;


protected:

    // Protected Member Functions

        //- Initialise containers and fields
        void initialise();

        //- Reset containers and fields
        void reset();


        // Evaluate

            //- Calculate force coefficients
            void calcForceCoeffs();

            //- Calculate moment coefficients
            void calcMomentCoeffs();

            //- Return integrated {total, pressure, viscous, porous} components
            List<scalar> integrateData(const List<Field<scalar>>& coeff) const;


        // Write

            //- Create the integrated-coefficient file
            bool createDataFile();

            //- Write header to the integrated-coefficient file
            bool writeDataFileHeader
            (
                const word& header,
                OFstream& os
            ) const;

            //- Write integrated coefficients to the integrated-coefficient file
            bool writeDataFile();

            //- Write integrated data to stream
            void logData
            (
                const word& descriptor,
                const List<scalar>& coeff
            ) const;


public:

    //- Runtime type information
    TypeName("forceCoeffs");


    // Constructors

        //- Construct from Time and dictionary
        forceCoeffs
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool readFields = true
        );


    //- Destructor
    virtual ~forceCoeffs() = default;


    // Member Functions

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
