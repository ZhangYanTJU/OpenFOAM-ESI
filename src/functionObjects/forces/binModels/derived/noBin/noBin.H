/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::binModels::noBin

Description
    Bin model providing an interface to bin properties without any active bin.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    <forces1/forceCoeffs1>
    {
        // Mandatory and optional entries
        ...

        // Optional entries (runtime modifiable)
        binModel    none;
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd | Deflt
      binModel     | Type name: noBin                   | word | no   | -
    \endtable

See also
  - Foam::binModel
  - Foam::binModels::singleDirectionUniformBin
  - Foam::binModels::uniformBin

SourceFiles
    noBin.C

\*---------------------------------------------------------------------------*/

#ifndef binModels_noBin_H
#define binModels_noBin_H

#include "binModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace binModels
{

/*---------------------------------------------------------------------------*\
                   Class noBin Declaration
\*---------------------------------------------------------------------------*/

class noBin
:
    public binModel
{
    // Generated Methods

        //- No copy construct
        noBin(const noBin&) = delete;

        //- No copy assignment
        void operator=(const noBin&) = delete;


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        noBin
        (
            const word& name,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~noBin() = default;


    // Member Functions

        // Access

            //- Read the dictionary
            virtual bool read(const dictionary& dict);


        // Evaluate

            //- Initialise bin properties
            virtual void initialise();

            //- Apply bins to force and moment containers
            virtual void applyBins
            (
                List<Field<vector>>& force,
                List<Field<vector>>& moment,
                const vectorField& d,
                const vectorField& Md,
                const vectorField& fN,
                const vectorField& fT
            );

            //- Apply bins to force and moment
            //- containers, and include porosity effects
            virtual void applyBins
            (
                List<Field<vector>>& force,
                List<Field<vector>>& moment,
                const vectorField& d,
                const vectorField& Md,
                const vectorField& fN,
                const vectorField& fT,
                const vectorField& fP
            );


        // Write

            //- Create binned-data files - no-op
            virtual void createBinnedDataFiles()
            {
                return;
            };

            //- Create a binned-data file - coefficients - no-op
            virtual void createBinnedDataFile
            (
                const word& title,
                const label id
            )
            {
                return;
            };

            //- Write binned force and moment data to files - no-op
            virtual void writeBinnedDataFiles
            (
                const List<Field<vector>>& force,
                const List<Field<vector>>& moment,
                const coordSystem::cartesian& coordSys
            )
            {
                return;
            };

            //- Write binned force and moment coefficient data to a file - no-op
            virtual void writeBinnedDataFile
            (
                List<Field<scalar>> coeff,
                const label id
            )
            {
                return;
            };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace binModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
