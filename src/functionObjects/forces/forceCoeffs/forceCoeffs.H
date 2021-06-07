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
    Extends the \c forces functionObject by providing coefficients for:
    - drag, side and lift forces (Cd, Cs, and Cl)
    - roll, pitch and yaw moments (CmRoll, CmPitch, and CmYaw)
    - front and rear axle force contributions (C(f) and C(r)) wherein

    \verbatim
        Cd(f/r) = 0.5*Cd \pm CmRoll
        Cs(f/r) = 0.5*Cs \pm CmYaw
        Cl(f/r) = 0.5*Cl \pm CmPitch
    \endverbatim

    The data can optionally be output into bins, defined in a given direction.

    The binned data provides the total and consitituent components per bin:
    - total coefficient
    - pressure coefficient contribution
    - viscous coefficient contribution
    - porous coefficient contribution

    Data is written into multiple files in the
    postProcessing/\<functionObjectName\> directory:
    - coefficient.dat   : integrated coefficients over all geometries
    - CdBin.dat         : drag coefficient bins
    - CsBin.dat         : side coefficient bins
    - ClBin.dat         : lift coefficient bins
    - CmRollBin.dat     : roll moment coefficient bins
    - CmPitchBin.dat    : pitch moment coefficient bins
    - CmYawBin.dat      : yaw moment coefficient bins

Usage
    Example of function object specification:
    \verbatim
    forceCoeffs1
    {
        type        forceCoeffs;
        libs        (forces);
        ...
        log         yes;
        writeFields yes;
        patches     (walls);

        // input keywords for directions of force/moment coefficients
        // refer below for options, and relations

        magUInf     100;
        lRef        3.5;
        Aref        2.2;
        porosity    no;

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
        Property     | Description                          | Required | Default
        type         | Type name: forceCoeffs               | yes |
        log          | Write force data to standard output  | no  | no
        writeFields  | Write force,moment coefficient fields | no | no
        patches      | Patches included in the forces calculation | yes |
        magUInf      | Free stream velocity magnitude       | yes |
        rhoInf       | Free stream density | for compressible cases |
        lRef         | Reference length scale for moment calculations | yes |
        Aref         | Reference area                       | yes |
        porosity     | Include porosity contributions       | no  | false
    \endtable

    Bin data is optional, but if the dictionary is present, the entries must
    be defined according to following:
    \table
        nBin         | Number of data bins                    | yes |
        direction    | Direction along which bins are defined | yes |
        cumulative   | Bin data accumulated with incresing distance | yes |
    \endtable

    Input of force/moment coefficient directions:
    - require an origin, and two orthogonal directions; the remaining orthogonal
    direction is determined accordingly.
    - can be added by the three options below.

    \verbatim
        CofR        (0 0 0); // Centre of rotation
        dragDir     (1 0 0);
        liftDir     (0 0 1);
    \endverbatim

    \verbatim
        origin (0 0 0);
        e1     (1 0 0);
        e3     (0 0 1); // combinations: (e1, e2) or (e2, e3) or (e3, e1)
    \endverbatim

    \verbatim
        coordinateSystem
        {
            origin  (0 0 0);
            rotation
            {
                type axes;
                e1 (1 0 0);
                e3 (0 0 1); // combinations: (e1, e2) or (e2, e3) or (e3, e1)
            }
        }
    \endverbatim

    The default direction relations are shown below:

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
    Foam::functionObject
    Foam::functionObjects::timeControl
    Foam::functionObjects::forces

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
