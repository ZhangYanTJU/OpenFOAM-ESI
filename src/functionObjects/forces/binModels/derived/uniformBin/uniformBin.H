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
    Foam::binModels::uniformBin

Description
    Bin model to divide a given patch/region into segments
    in given directions of a given Cartesian or cylindrical
    coordinate system, so that spatially-localised information
    (e.g. \c forces) can be output for each segment.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    <forces1/forceCoeffs1>
    {
        // Mandatory and ptional entries
        ...

        // Optional entries (runtime modifiable)
        binModel          uniformBin;

        binData
        {
            // Mandatory entries (runtime modifiable)
            nBin          (20 1 1);
            coordinateSystem
            {
                type      cylindrical;
                ...
            }

            // Optional entries (runtime modifiable)
            cumulative    false;
            minMax
            {
                e1    (5.5    1e300);
                e2    (-1e300 8);
                e3    (-1e300 1e300);
            }
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type | Reqd  | Deflt
      binModel  | Type name: uniformBin             | word | no    | -
      binData   | Bin settings             | dictionary | cndtnl       | -
      nBin      | Numbers of bins in (e1, e2, e3) directions | vector | yes | -
      coordinateSystem | Coordinate system settings | dictionary      | yes | -
      cumulative | Flag to bin data accumulated with increasing distance <!--
                 --> in a given binning direction   | bool   | no   | false
      minMax     | Min-max bounds in binning directions with respect to  <!--
                 --> the coordinateSystem's origin  | dictionary    | no    | -
    \endtable

    See the following for input coordinate-system examples:
     - \link coordinateSystem.H \endlink

Note
  - The order of bin numbering is (e1, e2, e3), where the first
    varies the fastest. For example, for a cylindrical bin with
    \f$ nBin = (radial, azimuth, height) = (2 4 2) \f$ the bin indices
    may look like as follows - note the counterclockwise increments:
    \verbatim
                           |-------------------|
                           | 12 |         | 14 |
                                | 11 | 13 |
                                | 9  | 15 |
                           | 10 |         | 16 |
                           |-------------------|
                          /    /         /    /
                         /    /         /    /
                        |-------------------|
                        | 4  |         | 6  |
                             | 3  | 5  |
                             | 1  | 7  |
                        | 2  |         | 8  |
                        |-------------------|
    \endverbatim

See also
  - Foam::binModel
  - Foam::binModels::noBin
  - Foam::binModels::singleDirectionUniformBin

SourceFiles
    uniformBin.C

\*---------------------------------------------------------------------------*/

#ifndef binModels_uniformBin_H
#define binModels_uniformBin_H

#include "binModel.H"
#include "writeFile.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace binModels
{

/*---------------------------------------------------------------------------*\
                   Class uniformBin Declaration
\*---------------------------------------------------------------------------*/

class uniformBin
:
    public binModel,
    public Foam::functionObjects::writeFile
{
protected:

        //- Flag to accumulate bin data with increasing
        //- distance in a given binning direction
        bool cumulative_;

        //- Numbers of bins in binning directions
        Vector<label> nBins_;

        //- Equidistant bin widths in binning directions
        vector binW_;

        //- Min-max bounds of bins in binning directions
        Vector<vector2D> binMinMax_;

        //- Local coordinate system of bins
        autoPtr<coordinateSystem> csysPtr_;


    // Generated Methods

        //- No copy construct
        uniformBin(const uniformBin&) = delete;

        //- No copy assignment
        void operator=(const uniformBin&) = delete;


public:

    //- Runtime type information
    TypeName("uniformBin");


    // Constructors

        //- Construct from components
        uniformBin
        (
            const word& name,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~uniformBin() = default;


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

            //- Create binned-data files
            virtual void createBinnedDataFiles();

            //- Create a binned-data file - coefficients
            virtual void createBinnedDataFile
            (
                const word& title,
                const label id
            );

            //- Write header for an binned-data file
            void writeBinnedDataFileHeader
            (
                const word& header,
                OFstream& os,
                const bool coeffStyle = false
            ) const;

            //- Write binned force and moment data to files
            virtual void writeBinnedDataFiles
            (
                const List<Field<vector>>& force,
                const List<Field<vector>>& moment,
                const coordSystem::cartesian& coordSys
            );

            //- Write binned force or moment data to a file
            void writeBinnedDataFile
            (
                List<Field<vector>>& fm,
                autoPtr<OFstream>& osPtr
            ) const;

            //- Write binned force and moment coefficient data to a file
            virtual void writeBinnedDataFile
            (
                List<Field<scalar>> coeff,
                const label id
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace binModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
