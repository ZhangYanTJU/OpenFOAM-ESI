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

#include "normBasicXiSubG.H"
#include "zeroGradientFvPatchField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace XiGModels
{
    defineTypeNameAndDebug(normBasicSubGrid, 0);
    addToRunTimeSelectionTable(XiGModel, normBasicSubGrid, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiGModels::normBasicSubGrid::normBasicSubGrid
(
    const dictionary& XiGProperties,
    const word& modelType,
    const psiuReactionThermo& thermo,
    const compressible::RASModel& turbulence,
    const volScalarField& Su
)
:
    XiGModel(XiGProperties, modelType, thermo, turbulence, Su),
//     Bv_
//     (
//         IOobject
//         (
//             "Bv",
//             Su.mesh().facesInstance(),
//             Su.mesh(),
//             IOobject::MUST_READ,
//             IOobject::NO_WRITE
//         ),
//         Su.mesh()
//     ),
    k1_(XiGModelCoeffs_.get<scalar>("k1")),
    kb1_(XiGModelCoeffs_.get<scalar>("kb1")),
    kbe_(XiGModelCoeffs_.get<scalar>("kbe")),
    kbx_(XiGModelCoeffs_.get<scalar>("kbx")),
    k2_(XiGModelCoeffs_.get<scalar>("k2")),
    LOverCw_(XiGModelCoeffs_.get<scalar>("LOverCw"))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::XiGModels::normBasicSubGrid::G() const
{
    const objectRegistry& db = Su_.db();
    const fvMesh& mesh = Su_.mesh();
    const volVectorField& U = db.lookupObject<volVectorField>("U");
    const volScalarField& b = db.lookupObject<volScalarField>("b");
    const volScalarField& Nv = db.lookupObject<volScalarField>("Nv");
    const volScalarField& St = db.lookupObject<volScalarField>("St");
    const volSymmTensorField& nsv = db.lookupObject<volSymmTensorField>("nsv");
    const volScalarField& Lobs = db.lookupObject<volScalarField>("Lobs");
    const volSymmTensorField& Bv = db.lookupObject<volSymmTensorField>("Bv");

    const scalarField Cw(pow(Su_.mesh().V(), 2.0/3.0));
    volScalarField CwVol
    (
        IOobject
        (
            "CwVol",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionSet(dimLength),
        Cw,
        zeroGradientFvPatchField<scalar>::typeName
    );
    CwVol.correctBoundaryConditions();

    if (!db.foundObject<volScalarField>("Ep"))
    {
        FatalErrorInFunction
            << "Looking for Ep in db that does not exist" << nl
            << Foam::abort(FatalError);
    }

    const volScalarField& Ep = db.lookupObject<volScalarField>("Ep");
    const volScalarField& Xp = db.lookupObject<volScalarField>("Xp");
    const volScalarField& Xi = db.lookupObject<volScalarField>("Xi");

    //tmp<volScalarField> tGtot = XiGModel_->G();
    auto tGtot = tmp<volScalarField>::New
    (
        IOobject
        (
            "tGtot",
            Su_.mesh().time().timeName(),
            Su_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        Su_.mesh(),
        dimensionedScalar(inv(dimTime), Zero)
    );
    auto& Gtot = tGtot.ref();

    //calculating flame normal

    volVectorField flNormal(fvc::reconstruct(fvc::snGrad(b)*mesh.magSf()));
    volScalarField mgb("mgb", mag(flNormal));
    dimensionedScalar dMgb("dMgb", mgb.dimensions(), SMALL);
    {
        volScalarField bc(b*(1.0-b));

        dMgb += 1.0e-8*
            (bc*mgb)().weightedAverage(mesh.V())
           /(bc.weightedAverage(mesh.V()) + SMALL);
    }
    mgb += dMgb;
    flNormal /= mgb;


    auto tN = tmp<volScalarField>::New
    (
        IOobject
        (
            "tN",
            Su_.mesh().time().timeName(),
            Su_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        Su_.mesh(),
        dimensionedScalar(Nv.dimensions(), Zero)
    );
    auto& N = tN.ref();

    N.primitiveFieldRef() = Nv.primitiveField()*Cw;

    auto tns = tmp<volSymmTensorField>::New
    (
        IOobject
        (
            "tns",
            Su_.mesh().time().timeName(),
            Su_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Su_.mesh(),
        dimensionedSymmTensor(nsv.dimensions(), Zero)
    );
    auto& ns = tns.ref();

    ns.primitiveFieldRef() = nsv.primitiveField()*Cw;

    auto tB = tmp<volSymmTensorField>::New
    (
        IOobject
        (
            "tB",
            Su_.mesh().time().timeName(),
            Su_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Su_.mesh(),
        dimensionedSymmTensor(Bv.dimensions(), Zero)
    );
    auto& B = tB.ref();

    B.primitiveFieldRef() = Bv.primitiveField()*sqrt(Cw);

    volScalarField Np(max(N - (flNormal & ns & flNormal), scalar(1)));

    // B_ is Bv*sqrt(Cw)
    volScalarField bl("bl",(flNormal & B & flNormal)/sqrt(Np));
    bl.min(1.0);

    volScalarField flSpeed("flSpeed", ((U & flNormal) + St)*b/(b+SMALL)) ;

    volScalarField up("up", sqrt((2.0/3.0)*turbulence_.k()));

    const volScalarField Gtot1
    (
        "Gtot1",
        (
            k1_ + kb1_*min(pow(bl, kbe_), kbx_)
        )*mag(flSpeed)/(max(Lobs, LOverCw_*CwVol))
    );

    const volScalarField Gtot2("Gtot2", k2_*Ep*Su_*Xi/(Xp - 0.999));

    const volScalarField value(pos(N - 1.e-3)*neg(bl - 0.99));

    Gtot = value*Gtot1+(1.0 - value)*Gtot2;


    /// if (Xi.mesh().time().outputTime())
    /// {
    ///     Gtot.write();
    ///     bl.write();
    ///     Lobs.write();
    ///     flSpeed.write();
    ///     N.write();
    /// }

    return tGtot;
}


Foam::tmp<Foam::volScalarField> Foam::XiGModels::normBasicSubGrid::Db() const
{
    // Not used //
    const objectRegistry& db = Su_.db();
    const volScalarField& Xi = db.lookupObject<volScalarField>("Xi");
    const volScalarField& rho = db.lookupObject<volScalarField>("rho");
    const volScalarField& mgb = db.lookupObject<volScalarField>("mgb");
    const volScalarField& Lobs = db.lookupObject<volScalarField>("Lobs");
    const volScalarField& Db = db.lookupObject<volScalarField>("Db");

    //return  turbulence_.muEff()
    return Db + rho*Su_*(Xi - 1.0)*mgb*(0.5*Lobs)*Lobs/(mgb*Lobs + 1.0);
}


bool Foam::XiGModels::normBasicSubGrid::read(const dictionary& XiGProperties)
{
    XiGModel::read(XiGProperties);

    XiGModelCoeffs_.readEntry("k1", k1_);
    XiGModelCoeffs_.readEntry("kb1", kb1_);
    XiGModelCoeffs_.readEntry("kbe", kbe_);
    XiGModelCoeffs_.readEntry("kbx", kbx_);
    XiGModelCoeffs_.readEntry("k2", k2_);

    return true;
}


// ************************************************************************* //
