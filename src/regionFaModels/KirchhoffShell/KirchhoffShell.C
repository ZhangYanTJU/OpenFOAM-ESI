/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "KirchhoffShell.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFields.H"
#include "zeroGradientFaPatchFields.H"
#include "subCycle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(KirchhoffShell, 0);

addToRunTimeSelectionTable(vibrationShellModel, KirchhoffShell, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool KirchhoffShell::read(const dictionary& dict)
{
    this->solution().readEntry("nNonOrthCorr", nNonOrthCorr_);
    return true;
}


void KirchhoffShell::solveDisplacement()
{
    if (debug)
    {
        InfoInFunction << endl;
    }
    
    const Time& time = primaryMesh().time();
    
    areaScalarField solidMass(rho()*h_);
    areaScalarField solidD(D()/solidMass);
    
    // Save old times
    areaScalarField w0 = w_.oldTime();
    areaScalarField w00 = w_.oldTime().oldTime();
    
    areaScalarField lap0(laplaceW_.oldTime());
    areaScalarField lapLap0(laplace2W_.oldTime());

    for
    (
        subCycleTime wSubCycle
        (
            const_cast<Time&>(time), 
            nSubCycles_
        );
        !(++wSubCycle).end();
    )
    {    
        laplaceW_ = fac::laplacian(w_);
        laplaceW_.correctBoundaryConditions();

        laplace2W_ = fac::laplacian(laplaceW_); 

        faScalarMatrix wEqn
        (
            fam::d2dt2(w_) 
          + f1_*fam::ddt(w_) 
          - f0_*sqrt(solidD)*fac::ddt(laplaceW_)
          + solidD*(laplace2W_ + f2_*fac::ddt(laplace2W_))
        ==
            ps_/solidMass
          + faOptions()(solidMass, w_, dimLength/sqr(dimTime))
        );
        
        faOptions().constrain(wEqn);

        wEqn.solve();
        w_.correctBoundaryConditions();
        w_.relax();
    }
    
    //wSubCycle.endSubCycle();


    // Correct index for w
    w_.timeIndex() = time.timeIndex();

    // Reset the old-time field value
    w_.oldTime() = w0;
    w_.oldTime().timeIndex() = time.timeIndex();
    w_.oldTime().oldTime() = w00;
    w_.oldTime().oldTime().timeIndex() = time.timeIndex();
    
    laplaceW_.oldTime() = lap0;
    laplaceW_.oldTime().timeIndex() = time.timeIndex();
    
    laplace2W_.oldTime() = lapLap0;
    laplace2W_.oldTime().timeIndex() = time.timeIndex();
    
    faOptions().correct(w_);
    
    Info<< "w min/max   = " << min(w_) << ", " << max(w_) << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

KirchhoffShell::KirchhoffShell
(
    const word& modelType,
    const fvPatch& patch,
    const dictionary& dict
)
:
    vibrationShellModel(modelType, patch, dict),
    f0_("f0", dimless, dict),
    f1_("f1", inv(dimTime), dict),
    f2_("f2", dimTime, dict),
    nNonOrthCorr_(1),
    nSubCycles_(1),
    ps_
    (
        IOobject
        (
            "ps_" + regionName_,
            primaryMesh().time().timeName(),
            primaryMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimPressure, Zero)
    ),
    h_
    (
        IOobject
        (
            "h_" + regionName_,
            primaryMesh().time().timeName(),
            primaryMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),
    laplaceW_
    (
        IOobject
        (
            "laplaceW_" + regionName_,
            primaryMesh().time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(inv(dimLength), Zero),
        w_.boundaryField().types()
    ),
    laplace2W_
    (
        IOobject
        (
            "laplace2W_" + regionName_,
            primaryMesh().time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(inv(pow3(dimLength)), Zero)
     )
{
    init();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

KirchhoffShell::~KirchhoffShell()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void KirchhoffShell::init()
{
}


void KirchhoffShell::preEvolveRegion()
{}


void KirchhoffShell::evolveRegion()
{
    nNonOrthCorr_ = solution().get<label>("nNonOrthCorr");
    nSubCycles_ = solution().get<label>("nSubCycles");

    for (int nonOrth=0; nonOrth<=nNonOrthCorr_; nonOrth++)
    {
        solveDisplacement();
    }
    
    // Update shell acceleration
    a_ = fac::d2dt2(w_);

   
}


const tmp<areaScalarField> KirchhoffShell::D() const
{
    const dimensionedScalar E("E", dimForce/dimArea , solid().E());
    const dimensionedScalar nu("nu", dimless, solid().nu());
    
    return tmp<areaScalarField>(E*pow3(h_)/(12*(1 - sqr(nu))));
}


const tmp<areaScalarField> KirchhoffShell::rho() const
{
    return tmp<areaScalarField>
    (
        new areaScalarField
        (
            IOobject
            (
                "rhos",
                primaryMesh().time().timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            regionMesh(),
            dimensionedScalar("rho", dimDensity, solid().rho()),
            zeroGradientFaPatchScalarField::typeName
        )
    );
}

void KirchhoffShell::info()
{
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
