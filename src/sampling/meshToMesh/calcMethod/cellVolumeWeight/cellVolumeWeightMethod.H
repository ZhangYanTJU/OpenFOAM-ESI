/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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
    Foam::cellVolumeWeightMethod

Description
    Cell-volume-weighted mesh-to-mesh interpolation class

    Volume conservative.

SourceFiles
    cellVolumeWeightMethod.C

\*---------------------------------------------------------------------------*/

#ifndef cellVolumeWeightMethod_H
#define cellVolumeWeightMethod_H

#include "meshToMeshMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class cellVolumeWeightMethod Declaration
\*---------------------------------------------------------------------------*/

class cellVolumeWeightMethod
:
    public meshToMeshMethod
{
protected:

    // Protected Member Functions

        //- Find indices of overlapping cells in src and tgt meshes - returns
        //  true if found a matching pair
        bool findInitialSeeds
        (
            const labelList& srcCellIDs,
            const boolList& mapFlag,
            const label startSeedI,
            label& srcSeedI,
            label& tgtSeedI
        ) const;

        //- Calculate the mesh-to-mesh addressing and weights
        void calculateAddressing
        (
            labelListList& srcToTgtCellAddr,
            scalarListList& srcToTgtCellWght,
            labelListList& tgtToSrcCellAddr,
            scalarListList& tgtToSrcCellWght,
            const label srcSeedI,
            const label tgtSeedI,
            const labelList& srcCellIDs,
            boolList& mapFlag,
            label& startSeedI
        );

        //- Set the next cells in the advancing front algorithm
        void setNextCells
        (
            label& startSeedI,
            label& srcCelli,
            label& tgtCelli,
            const labelList& srcCellIDs,
            const boolList& mapFlag,
            const labelUList& visitedCells,
            labelList& seedCells
        ) const;


private:

    // Private member functions

        //- No copy construct
        cellVolumeWeightMethod(const cellVolumeWeightMethod&) = delete;

        //- No copy assignment
        void operator=(const cellVolumeWeightMethod&) = delete;


public:

    //- Run-time type information
    TypeName("cellVolumeWeight");

    //- Construct from source and target meshes
    cellVolumeWeightMethod(const polyMesh& src, const polyMesh& tgt);

    //- Destructor
    virtual ~cellVolumeWeightMethod();


    // Member Functions

        // Evaluate

            //- Calculate addressing and weights and optionally offset vectors
            virtual void calculate
            (
                labelListList& srcToTgtAddr,
                scalarListList& srcToTgtWght,
                pointListList& srcToTgtVec,
                labelListList& tgtToSrcAddr,
                scalarListList& tgtToSrcWght,
                pointListList& tgtToSrcVec
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
