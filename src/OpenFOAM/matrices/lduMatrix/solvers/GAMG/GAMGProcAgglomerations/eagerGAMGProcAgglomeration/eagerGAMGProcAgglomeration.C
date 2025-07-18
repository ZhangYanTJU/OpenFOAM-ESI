/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "eagerGAMGProcAgglomeration.H"
#include "addToRunTimeSelectionTable.H"
#include "GAMGAgglomeration.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(eagerGAMGProcAgglomeration, 0);

    addToRunTimeSelectionTable
    (
        GAMGProcAgglomeration,
        eagerGAMGProcAgglomeration,
        GAMGAgglomeration
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eagerGAMGProcAgglomeration::eagerGAMGProcAgglomeration
(
    GAMGAgglomeration& agglom,
    const dictionary& controlDict
)
:
    GAMGProcAgglomeration(agglom, controlDict),
    mergeLevels_(controlDict.getOrDefault<label>("mergeLevels", 1))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eagerGAMGProcAgglomeration::~eagerGAMGProcAgglomeration()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::eagerGAMGProcAgglomeration::agglomerate()
{
    if (debug)
    {
        Pout<< nl << "Starting mesh overview" << endl;
        printStats(Pout, agglom_);
    }

    if (agglom_.size() >= 1)
    {
        // Agglomerate one but last level (since also agglomerating
        // restrictAddressing)
        for
        (
            label fineLevelIndex = 2;
            fineLevelIndex < agglom_.size();
            fineLevelIndex++
        )
        {
            if (agglom_.hasMeshLevel(fineLevelIndex))
            {
                // Get the fine mesh
                const lduMesh& levelMesh = agglom_.meshLevel(fineLevelIndex);
                label levelComm = levelMesh.comm();
                label nProcs = UPstream::nProcs(levelComm);

                if (nProcs > 1)
                {
                    // Processor restriction map: per processor the coarse
                    // processor
                    labelList procAgglomMap(nProcs);

                    forAll(procAgglomMap, proci)
                    {
                        procAgglomMap[proci] = proci/(1<<mergeLevels_);
                    }

                    // Master processor
                    labelList masterProcs;
                    // Local processors that agglomerate. agglomProcIDs[0]
                    // is in masterProc.
                    List<label> agglomProcIDs;
                    GAMGAgglomeration::calculateRegionMaster
                    (
                        levelComm,
                        procAgglomMap,
                        masterProcs,
                        agglomProcIDs
                    );

                    // Communicator for the processor-agglomerated matrix
                    comms_.push_back
                    (
                        UPstream::newCommunicator
                        (
                            levelComm,
                            masterProcs
                        )
                    );

                    // Use processor agglomeration maps to do the actual
                    // collecting.
                    if (UPstream::myProcNo(levelComm) != -1)
                    {
                        GAMGProcAgglomeration::agglomerate
                        (
                            fineLevelIndex,
                            procAgglomMap,
                            masterProcs,
                            agglomProcIDs,
                            comms_.back()
                        );
                    }
                }
            }
        }
    }

    // Print a bit
    if (debug)
    {
        Pout<< nl << "Agglomerated mesh overview" << endl;
        printStats(Pout, agglom_);
    }

    return true;
}


// ************************************************************************* //
