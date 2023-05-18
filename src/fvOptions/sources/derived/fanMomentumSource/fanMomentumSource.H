/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Louis Vittoz, SimScale GmbH
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
    Foam::fanMomentumSource

Description
    This source term models the action of a fan on the flow.

Usage
    \table
        Property     | Description                          | Required | Default
        fanCurve     | Pressure drop vs flow-rate           | yes |
        fields       | Name of operand field                | yes |
        flowDir      | Direction of the flow through the fan| yes |
        faceZone     | Name of upstream faceZone            | yes |
        thickness    | Thickness of the fan                 | yes |
     \endtable

    The inherited entries are elaborated in:
     - \link cellSetOption.H \endlink

    Example of the source term specification:
    \verbatim
    fanModel
    {
        type          fanMomentumSource;
        active        on;
        selectionMode cellZone;
        cellZone      fanCellZone;
        faceZone      fanUpstreamFaceZone;
        fanCurve
        {
            type table;
            file "constant/fanCurve";
        }
        fields   (U);
        flowDir  (1.0 0.0 0.0);
        thickness 0.2;
    }

Note
    The fan curve needs to provide for a pressure drop expressed in Pascal and
    is specified as a function of the volumetric flow rate in (m^3/s)

SourceFiles
    fanMomentumSource.C


\*---------------------------------------------------------------------------*/

#ifndef fanMomentumSource_H
#define fanMomentumSource_H

#include "autoPtr.H"
#include "fvMesh.H"
#include "volFields.H"
#include "cellSetOption.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class fanMomentumSource Declaration
\*---------------------------------------------------------------------------*/

class fanMomentumSource
:
    public fv::cellSetOption
{
    // Private Data

        //- Volumetric flow rate vs. pressure drop table
        autoPtr<Function1<scalar>> fanCurve_;

        //- Direction of flow through the fan
        vector flowDir_;

        //- Thickness of the fan
        const scalar thickness_;

        //- Pressure drop
        scalar gradPFan_;

        //- Local list of upstream face IDs
        labelList upstreamFaceIDs_;

        //- Local list of upstream patch ID per face
        labelList upstreamFacePatchIDs_;

        //- Id for the surrounding face zone
        label surroundingFaceZoneID_;

        //- Cells in zone
        labelHashSet cellsInZones_;

        //- Reference density for incompressible cases
        autoPtr<scalar> rho_;


    // Private Member Functions

        //- Write the pressure gradient to file (for restarts etc)
        void writeProps(const scalar gradP, const scalar flowRate) const;

        //- Collect all faces upstream of the centre of gravity of the cell
        //  zone
        void initializeUpstreamFaces();

        //- Calculate the volumetric flow rate
        template<typename flowRateFunctorPatch, typename flowRateFunctor>
        scalar calculateFlowRate
        (
            flowRateFunctorPatch fPatch,
            flowRateFunctor f
        ) const;


public:

    //- Runtime type information
    TypeName("fanMomentumSource");


    // Constructors

        //- Construct from explicit source name and mesh
        fanMomentumSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        fanMomentumSource(const fanMomentumSource&) = delete;

        //- No copy assignment
        void operator=(const fanMomentumSource&) = delete;


    //- Destructor
    virtual ~fanMomentumSource() = default;


    // Member Functions

        // Evaluation

            //- Add explicit contribution to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );

        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
