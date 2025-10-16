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

Class
    Foam::functionObjects::radiometerProbes

Group
    grpUtilitiesFunctionObjects

Description
    Probes the incident radiative heat flux, qin, at arbitrary points within a
    domain.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    radiometer
    {
        // Mandatory entries
        type                radiometerProbes;
        libs                (utilityFunctionObjects);
        probeLocations      (<vectorList>);
        probeNormals        (<vectorList>);

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd | Deflt
      type         | Type name: radiometerProbes        | word | yes  | -
      libs         | Library name: utilityFunctionObjects | word | yes  | -
      probeLocations | Locations of probes              | vectorList | yes | -
      probeNormals | Normals of specified probes        | vectorList | yes | -
    \endtable

    The inherited entries are elaborated in:
     - \link regionFunctionObject.H \endlink
     - \link internalProber.H \endlink
     - \link writeFile.H \endlink
     - \link fvDOM.H \endlink

Note
  - The function object can only be used with the \c fvDOM radiation model.
  - The \c solverFreq input of the \c fvDOM model has superiority over
    \c executeControl and \c writeControl entries.

SourceFiles
    radiometerProbes.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_radiometerProbes_H
#define Foam_functionObjects_radiometerProbes_H

#include "regionFunctionObject.H"
#include "internalProber.H"
#include "writeFile.H"
#include "fvDOM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                     Class radiometerProbes Declaration
\*---------------------------------------------------------------------------*/

class radiometerProbes
:
    public regionFunctionObject,
    public internalProber,
    public writeFile
{
    // Private Data

        //- Const reference to the underlying radiation model
        const radiation::fvDOM& dom_;

        //- Normal vectors of the specified probes
        vectorField n_;

        //- Pre-computed inner product of probe normals (n_) and average
        //- solid-angle direction (dAve) per radiative intensity ray
        List<scalarList> n_dAve_;

        //- Directional selection coefficient for radiative intensity rays
        //  false: ray entering the probe
        //  true: ray leaving the probe
        List<boolList> C_;

        //- Incident radiative heat flux per probe location
        scalarField qin_;

        //- Number of radiative intensity rays
        label nRay_;

        //- Number of probe locations/normals
        label szProbes_;

        //- Flag to identify whether the iteration is the first iteration
        //  Resets with a restarted simulation
        bool firstIter_;


    // Private Member Functions

        //- Write file-header information into the output file
        virtual void writeFileHeader(Ostream& os);

        //- Return the flag to decide if radiation-model calculations are
        //- performed, so that function object calculations can proceed
        bool shouldCalcThisStep() const
        {
            return
                firstIter_
             || (mesh_.time().timeIndex() % dom_.solverFreq() == 0);
        }


public:

    //- Runtime type information
    TypeName("radiometerProbes");


    // Generated Methods

        //- No copy construct
        radiometerProbes(const radiometerProbes&) = delete;

        //- No copy assignment
        void operator=(const radiometerProbes&) = delete;


    // Constructors

        //- Construct from name, Time and dictionary
        radiometerProbes
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~radiometerProbes() = default;


    // Public Member Functions

        //- Read the function object settings
        virtual bool read(const dictionary&);

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
