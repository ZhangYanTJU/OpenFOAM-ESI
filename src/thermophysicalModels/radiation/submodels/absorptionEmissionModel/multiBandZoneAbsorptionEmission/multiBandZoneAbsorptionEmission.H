/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2019 OpenCFD Ltd.
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
    Foam::radiation::multiBandZoneAbsorptionEmission

Group
    grpRadiationAbsorptionEmissionSubModels

Description

    multiBandZoneAbsorptionEmission radiation absorption/emission for solids.


SourceFiles
    multiBandZoneAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef multiBandZoneAbsorptionEmission_H
#define multiBandZoneAbsorptionEmission_H

#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class multiBandZoneAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class multiBandZoneAbsorptionEmission
:
    public absorptionEmissionModel
{
public:

    // Public data

        //- Maximum number of bands
        static const int maxBands_ = 5;


private:

    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //-  Absorption coefficients outside cellSets
        scalarList absCoeffs_;

        //-  Emissivity coefficients outside cellSets
        scalarList emiCoeffs_;

        //- Bands
        label nBands_;

        //- Cell zones absorptivity
        HashTable<scalarList> zoneAbsorptivity_;

        //- Cell zones emissivity
        HashTable<scalarList> zoneEmissivity_;

        //- The cellZones selected
        labelList zoneIds_;



public:

    //- Runtime type information
    TypeName("multiBandZoneAbsorptionEmission");


    // Constructors

        //- Construct from components
        multiBandZoneAbsorptionEmission
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~multiBandZoneAbsorptionEmission();


    // Member Functions

        // Access

            // Absorption coefficient

                //- Absorption coefficient
                tmp<volScalarField> aCont(const label bandI) const;


            // Emission coefficient

                //- Emission coefficient
                tmp<volScalarField> eCont(const label bandI) const;


            // Emission contribution

                //- Emission contribution
                tmp<volScalarField> ECont(const label bandI) const;


        inline bool isGrey() const
        {
            return false;
        }

        //- Number of bands
        inline label nBands() const
        {
            return nBands_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
