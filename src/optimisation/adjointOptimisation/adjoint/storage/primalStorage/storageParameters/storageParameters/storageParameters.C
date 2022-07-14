/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 PCOpt/NTUA
    Copyright (C) 2022      FOSS GP
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

#include "Ostream.H"
#include "storageParameters.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(storageParameters, 0);
defineRunTimeSelectionTable(storageParameters, dictionary);


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

label storageParameters::counter = 0;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void storageParameters::findCellsNextToWall()
{
    // Find the boundary patches that are of type wall
    const fvPatchList& patches = mesh_.boundary();
    DynamicList<label> patchList;
    label size(0);
    forAll(patches, patchI)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        if (isA<wallFvPatch>(patch))
        {
            patchList.append(patchI);
            size += patch.faceCells().size();
        }
    }
    // Build list containing the IDs of the first cell centers next to wall
    // boundaries
    wallList_.setSize(size, -1);
    label pos(0);
    forAll(patchList, patchI)
    {
        const fvPatch& patch = mesh_.boundary()[patchList[patchI]];
        forAll(patch, faceI)
        {
            wallList_[pos++] = patch.faceCells()[faceI];
        }
    }
    Info<< "Found " << patchList.size() << " wall boundary(ies) "
        << "with total " << returnReduce(wallList_.size(), sumOp<label>())
        << " cells next to them." << endl;
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<storageParameters> storageParameters::New
(
    const fvMesh& mesh,
    const dictionary& storageDict,
    const wordList& allocatedFieldNames
)
{
    const dictionary& dict = storageDict.subDict("storageParameters");
    const word type(dict.get<word>("algorithm") + "Parameters" );
    Info<< "Compression parameters type" << tab << type << endl;

    auto* ctorPtr = dictionaryConstructorTable(type);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            storageDict,
            "storageParameters",
            type,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<storageParameters>
        (ctorPtr(mesh, storageDict, allocatedFieldNames));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

storageParameters::storageParameters
(
    const fvMesh& mesh,
    const dictionary& storageDict,
    const wordList& allocatedFieldNames
)
:
    mesh_(mesh),
    dict_(storageDict.subDict("storageParameters")),
    allocatedFieldNames_(allocatedFieldNames),
    compressionMethod_(),
    writeFieldTimes_(),
    timing_(true),
    writeAll_(false),
    tiny_(1.e-18),
    solDirs_((mesh_.solutionD() + Vector<label>::one)/2)
{
    if (Pstream::master())
    {
        counter++;
    }
    initialize();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void storageParameters::initialize()
{
    // Build list containing the IDs of the first cell centers next to wall
    // boundaries
    findCellsNextToWall();

    adjustTimeStep_ =
        mesh_.time().controlDict().getOrDefault<Switch>("adjustTimeStep", false);
    storeAllBoundaries_ =
        dict_.getOrDefault<Switch>("storeAllBoundaries", false);
    storeUniformBoundaries_ =
        dict_.getOrDefault<Switch>("storeUniformBoundaries", false);
    algorithm_ = dict_.get<word>("algorithm");
    timing_ = dict_.getOrDefault<bool>("timing", true);
    Info<< nl << "General Parameters:" << nl
        << tab << "algorithm              " << algorithm_ << nl
        << tab << "timing                 " << timing_ << nl
        << tab << "storeAllBoundaries     " << storeAllBoundaries_ << endl;
    if (!timing_)
    {
        WarningInFunction
            << "Timing should be deactivated only when a reduced performance "
            << "is acceptable."
            << endl;
    }
    if (!storeAllBoundaries_)
    {
        Info<< tab << "storeUniformBoundaries "
            << storeUniformBoundaries_ << endl;
        if (storeUniformBoundaries_)
        {
            WarningInFunction
                << "'storeUniformBoundaries' has been activated." << nl
                << "This is redundant but useful in debugging." << nl
                << "Make sure it is intentional." << endl;
        }
    }
    if (storeAllBoundaries_)
    {
        WarningInFunction
            << "'storeAllBoundaries' should be activated only for debugging "
            << "purposes, since unneccessary storage is allocated" << endl;
    }
    writeAll_ = dict_.getOrDefault<bool>("writeAll", false);
    Info << tab << "writeAll               " << writeAll_ << endl;
    if (!writeAll_)
    {
        writeFieldTimes_ =
            dict_.getOrDefault<scalarList>("writeFieldTimes", scalarList());
        if (writeFieldTimes_.size() > 0)
        {
            Info<< tab << "Time-steps to write primal fields:" << endl
                << tab << writeFieldTimes_ << nl << endl;
        }
    }
    Info << endl;
}


void storageParameters::reset()
{
    // Does nothing in base
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
