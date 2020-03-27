/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 DLR
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

#include "reconstructionSchemes.H"
#include "OFstream.H"
#include "cutCellPLIC.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(reconstructionSchemes, 0);
    defineRunTimeSelectionTable(reconstructionSchemes, components);
}



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::reconstructionSchemes::alreadyReconstructed()
{
    const fvMesh& mesh = alpha1_.mesh();
    label& curTimeIndex = timeIndexAndIter_.first();
    label& curIter = timeIndexAndIter_.second();

    // rest timeIndex and curIter
    if (mesh.time().timeIndex() > curTimeIndex)
    {
        curTimeIndex = mesh.time().timeIndex();
        curIter = 0;
        return false;
    }

    // reconstruct always when subcycling
    if (mesh.time().subCycling() != 0)
    {
        return false;
    }

    ++curIter;
    if (curIter > 1)
    {
        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reconstructionSchemes::reconstructionSchemes
(
    const word& type,
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U,
    const dictionary& dict
)
:
    IOdictionary
    (
        IOobject
        (
            "reconstructionScheme",
            alpha1.time().constant(),
            alpha1.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    reconstructionSchemesCoeffs_(dict),
    alpha1_(alpha1),
    phi_(phi),
    U_(U),
    normal_
    (
        IOobject
        (
            "recon::normal",
            alpha1_.mesh().time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedVector(dimArea, Zero)
    ),
    centre_
    (
        IOobject
        (
            "recon::centre",
            alpha1_.mesh().time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedVector(dimLength, Zero)
    ),
    interfaceCell_(alpha1_.mesh().nCells(), false),
    interfaceLabels_(0.2*alpha1_.mesh().nCells()),
    timeIndexAndIter_(0, 0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dictionary& Foam::reconstructionSchemes::modelDict() const
{
    return reconstructionSchemesCoeffs_;
}


Foam::dictionary& Foam::reconstructionSchemes::modelDict()
{
    return reconstructionSchemesCoeffs_;
}

Foam::reconstructionSchemes::interface Foam::reconstructionSchemes::surface()
{
    reconstruct();
    const fvMesh& mesh = centre_.mesh();

    cutCellPLIC cellCut(mesh);

    DynamicList<List<point>> facePts;
    DynamicList<label> interfaceCellAdressing(0.1*mesh.nCells());

    forAll(interfaceCell_,cellI)
    {
        if (interfaceCell_[cellI])
        {
            if (mag(normal_[cellI]) != 0)
            {
                interfaceCellAdressing.append(cellI);
                vector n = -normal_[cellI]/mag(normal_[cellI]);

                scalar cutVal = (centre_[cellI]-mesh.C()[cellI]) & n;

                cellCut.calcSubCell(cellI,cutVal,n);
                facePts.append(cellCut.facePoints());
            }
        }
    }

    labelList meshCells(interfaceCellAdressing.size());

    forAll(meshCells,i)
    {
        meshCells[i] = interfaceCellAdressing[i];
    }

    // Transfer to mesh storage
    
    faceList faces(facePts.size());

    label nPoints = 0;
    forAll(facePts,i)
    {
        face f(facePts[i].size());
        forAll(f,fi)
        {
            f[fi] = nPoints + fi;
        }
        faces[i] = f;

        nPoints += facePts[i].size();
    }
    pointField pts(nPoints);

    nPoints = 0; // reuse
    forAll(facePts,i)
    {
        forAll(facePts[i],fi)
        {
            pts[nPoints] = facePts[i][fi];
            ++nPoints;
        }
    }

    interface surf(pts, faces, meshCells);

    return surf;
}


// ************************************************************************* //
