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

    DynamicList<point> dynPts;
    DynamicList<face> dynFaces(mesh.nCells()*0.1);
    bitSet interfaceCellAddressing(mesh.nCells());

    forAll(interfaceCell_,celli)
    {
        if (interfaceCell_[celli])
        {
            if (mag(normal_[celli]) != 0)
            {
                interfaceCellAddressing.set(celli);
                vector n = -normal_[celli]/mag(normal_[celli]);

                scalar cutVal = (centre_[celli]-mesh.C()[celli]) & n;
                cellCut.calcSubCell(celli,cutVal,n);
                // cellCut.facePoints() are ordered and not
                // connected to the other face 
                // append face with the starting label: dynPts.size()
                dynFaces.append
                (
                    face(identity(cellCut.facePoints().size(), dynPts.size()))
                );
                dynPts.append(cellCut.facePoints());
            }
        }
    }

    labelList meshCells(interfaceCellAddressing.toc());

    // Transfer to mesh storage
    pointField pts(std::move(dynPts));
    faceList faces(std::move(dynFaces));

    interface surf(pts, faces, meshCells);

    return surf;
}


// ************************************************************************* //
