/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
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
    Foam::binModels::singleDirectionUniformBin

Description
    Bin model to divide a given patch/region into segments
    in a given single direction, so that spatially-localised
    information (e.g. \c forces) can be output for each segment.

    For example, for a 10cm-long patch extending only in the x-direction
    can be binned into 5 bins in the x-direction, so that local information
    are output for each 2cm-long segment.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    <forces1/forceCoeffs1>
    {
        // Mandatory and other optional entries
        ...

        // Optional (inherited) entries (runtime modifiable)
        binModel          singleDirectionUniformBin;

        // Conditional optional entries (runtime modifiable)
        binData
        {
            // Mandatory entries (runtime modifiable)
            nBin          5;
            direction     (1 0 0);
            cumulative    true;

            // Optional entries (runtime modifiable)
            min           0.0;
            max           2.5;
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                           | Type | Reqd  | Dflt
      binModel  | Type name: singleDirectionUniformBin  | word | no    | -
      binData   | Bin settings             | dictionary | cndtnl       | -
      nBin      | Number of bins in the given direction | vector | yes | -
      direction | Direction used to determine bin orientation | vector | yes | -
      cumulative | Flag to bin data accumulated with increasing distance <!--
                 --> in a given binning direction       | bool   | yes | -
      min       | Min-bound in the binning direction with respect to <!--
                 --> the global coordinate system's origin | scalar | no | GREAT
      max       | Max-bound in the binning direction with respect to <!--
                 --> the global coordinate system's origin | scalar | no | GREAT
    \endtable

Note
  - \c singleDirectionUniformBin is the same bin method available
    in the pre-v2012 versions. To ensure backward compatibility,
    \c singleDirectionUniformBin bin model is automatically activated
    when only \c binData entry is present within the function-object
    dictionaries of \c forces and \c forceCoeffs.

See also
  - Foam::binModel
  - Foam::binModels::noBin
  - Foam::binModels::uniformBin

SourceFiles
    singleDirectionUniformBin.C

\*---------------------------------------------------------------------------*/

#ifndef binModels_singleDirectionUniformBin_H
#define binModels_singleDirectionUniformBin_H

#include "binModel.H"
#include "writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace binModels
{

/*---------------------------------------------------------------------------*\
                   Class singleDirectionUniformBin Declaration
\*---------------------------------------------------------------------------*/

class singleDirectionUniformBin
:
    public binModel
{
protected:

        //- Flag to accumulate bin data with increasing
        //- distance in a given binning direction
        bool cumulative_;

        //- Distance between bin divisions
        scalar binDx_;

        //- Minimum bin bound
        scalar binMin_;

        //- Maximum bin bound
        scalar binMax_;

        //- Direction used to determine bin orientation
        vector binDir_;


    // Protected Member Functions

        //- Write header for a binned-data file
        template<class Type>
        void writeFileHeader
        (
            OFstream& os,
            const bool coeffStyle = false
        ) const;

        //- Initialise bin properties
        virtual void initialise();

        //- Apply the binning to field fieldi
        template<class Type>
        bool processField(const label fieldi);

        //- Write binned data to ostream
        template<class Type>
        void writeBinnedData
        (
            List<List<Type>>& data,
            Ostream& os
        ) const;


public:

    //- Runtime type information
    TypeName("singleDirectionUniformBin");


    // Constructors

        //- Construct from components
        singleDirectionUniformBin
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& outputPrefix
        );

        //- No copy construct
        singleDirectionUniformBin(const singleDirectionUniformBin&) = delete;

        //- No copy assignment
        void operator=(const singleDirectionUniformBin&) = delete;


    //- Destructor
    virtual ~singleDirectionUniformBin() = default;


    // Member Functions

        //- Read the dictionary
        virtual bool read(const dictionary& dict);

        //- Apply bins to force and moment containers
        virtual void apply();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace binModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "singleDirectionUniformBinTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
