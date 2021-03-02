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

\*---------------------------------------------------------------------------*/

#include "normBasicXiSubXiEq.H"
#include "addToRunTimeSelectionTable.H"
#include "ignition.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace XiEqModels
{
    defineTypeNameAndDebug(normBasicSubGrid, 0);
    addToRunTimeSelectionTable(XiEqModel, normBasicSubGrid, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiEqModels::normBasicSubGrid::normBasicSubGrid
(
    const dictionary& XiEqProperties,
    const word& modelType,
    const psiuReactionThermo& thermo,
    const compressible::RASModel& turbulence,
    const volScalarField& Su
)
:
    XiEqModel(XiEqProperties, modelType,thermo, turbulence, Su),
    Cxpe1_(XiEqModelCoeffs_.get<scalar>("Cxpe1")),
    Cxpe2_(XiEqModelCoeffs_.get<scalar>("Cxpe2")),
    Cxpe3_(XiEqModelCoeffs_.get<scalar>("Cxpe3")),
    Cxpe4_(XiEqModelCoeffs_.get<scalar>("Cxpe4"))
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::XiEqModels::normBasicSubGrid::~normBasicSubGrid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::XiEqModels::normBasicSubGrid::XiEq() const
{
    const fvMesh& mesh = Su_.mesh();
    const volVectorField& U = mesh.lookupObject<volVectorField>("U");
    const volScalarField& b = mesh.lookupObject<volScalarField>("b");
    const volScalarField& Nv = mesh.lookupObject<volScalarField>("Nv");
    const volSymmTensorField& nsv =
        mesh.lookupObject<volSymmTensorField>("nsv");
    const volSymmTensorField& Bv =
        mesh.lookupObject<volSymmTensorField>("Bv");

    volScalarField magU(mag(U));

    const scalarField Cw = pow(mesh.V(), 2.0/3.0);

    tmp<volScalarField> tN
    (
        new volScalarField
        (
            IOobject
            (
                "tN",
                mesh.time().constant(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(Nv.dimensions(), Zero)
        )
    );

    volScalarField& N = tN.ref();

    N.primitiveFieldRef() = Nv.primitiveField()*Cw;

    tmp<volSymmTensorField> tns
    (
        new volSymmTensorField
        (
            IOobject
            (
                "tns",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedSymmTensor(nsv.dimensions(), Zero)
        )
    );
     volSymmTensorField& ns = tns.ref();

    tmp<volSymmTensorField> tB
    (
        new volSymmTensorField
        (
            IOobject
            (
                "tB",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedSymmTensor(Bv.dimensions(), Zero)
        )
    );
    volSymmTensorField& B = tB.ref();

    //calculating flame normal

    volVectorField flNormal
    (
        "flNormal",
        fvc::reconstruct(fvc::snGrad(b)*mesh.magSf())
    );

    volScalarField mgb("mgb", mag(flNormal));

    dimensionedScalar dMgb("dMgb", mgb.dimensions(), SMALL);

    const volScalarField bc(b*(1.0-b));

    dMgb += 1.0e-8*
        (bc*mgb)().weightedAverage(mesh.V())
       /(bc.weightedAverage(mesh.V()) + SMALL);

    mgb += dMgb;
    flNormal /= mgb;

    B.primitiveFieldRef() = Bv.primitiveField()*sqrt(Cw);
    volScalarField Ntemp("Ntemp", N);
    volScalarField Np("Np", max(N - (flNormal & ns & flNormal), scalar(1)));

    // B_ is Bv*sqrt(Cw)
    volScalarField bl("bl",(flNormal & B & flNormal)/sqrt(Np));
    bl.min(1.0);

    volScalarField up(sqrt((2.0/3.0)*turbulence_.k()));

    IOdictionary combustionProperties
    (
        IOobject
        (
            "combustionProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    ignition ign(combustionProperties, mesh.time(), mesh);

    dimensionedVector ignLoc("ignLoc", dimLength, ign.sites()[0].location());

    dimensionedScalar filtRad2
    (
        "filtRad2",
        dimLength,
        6.0*ign.sites()[0].diameter()
    );

    const volScalarField filDist(mag(mesh.C() - ignLoc));

    const volScalarField filterMult
    (
        pos(filDist - filtRad2)*neg(bl - 0.99)*pos(N - 1e-3)
    );

    tmp<volScalarField> XiSubEq
    (
        scalar(1)
      + min( min(Cxpe1_, Cxpe2_*magU/up)*sqrt(bl), Cxpe3_)
      * filterMult
    );

    return XiSubEq;
}


bool Foam::XiEqModels::normBasicSubGrid::read(const dictionary& XiEqProperties)
{
    XiEqModel::read(XiEqProperties);

    XiEqModelCoeffs_.readEntry("Cxpe1", Cxpe1_);
    XiEqModelCoeffs_.readEntry("Cxpe2", Cxpe2_);
    XiEqModelCoeffs_.readEntry("Cxpe3", Cxpe3_);
    XiEqModelCoeffs_.readEntry("Cxpe4", Cxpe4_);

    return true;
}


// ************************************************************************* //
