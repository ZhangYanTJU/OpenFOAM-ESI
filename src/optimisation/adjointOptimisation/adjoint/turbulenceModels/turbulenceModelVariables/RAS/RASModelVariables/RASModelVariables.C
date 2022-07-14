/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2022 PCOpt/NTUA
    Copyright (C) 2013-2022 FOSS GP
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "RASModelVariables.H"
#include "addToRunTimeSelectionTable.H"
#include "variablesSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(RASModelVariables, 0);
defineRunTimeSelectionTable(RASModelVariables, dictionary);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void RASModelVariables::allocateInitValues()
{
    if (solverControl_.storeInitValues())
    {
        Info<< "Storing initial values of turbulence variables" << endl;

        if (hasTMVar1())
        {
            variablesSet::setInitField(TMVar1InitPtr_, TMVar1Inst());
        }

        if (hasTMVar2())
        {
            variablesSet::setInitField(TMVar2InitPtr_, TMVar2Inst());
        }

        if (hasNut())
        {
            variablesSet::setInitField(nutInitPtr_, nutRefInst());
        }
    }
}


void RASModelVariables::allocateMeanFields()
{
    if (solverControl_.average())
    {
        Info<< "Allocating mean values of turbulence variables" << endl;

        if (hasTMVar1())
        {
            variablesSet::setMeanField(TMVar1MeanPtr_, TMVar1Inst(), mesh_);
        }

        if (hasTMVar2())
        {
            variablesSet::setMeanField(TMVar2MeanPtr_, TMVar2Inst(), mesh_);
        }

        if (hasNut())
        {
            variablesSet::setMeanField(nutMeanPtr_, nutRefInst(), mesh_);
        }
    }
}


Foam::refPtr<Foam::volScalarField>
RASModelVariables::cloneRefPtr(const refPtr<volScalarField>& obj) const
{
    if (obj)
    {
        const volScalarField& sf = obj();

        const word timeName = mesh_.time().timeName();

        return refPtr<volScalarField>::New(sf.name() + timeName, sf);
    }

    return nullptr;
}


void RASModelVariables::copyAndRename
(
    volScalarField& f1,
    volScalarField& f2
)
{
    f1 == f2;
    const word name1 = f1.name();
    const word name2 = f2.name();

    // Extra rename to avoid databese collision
    f2.rename("temp");
    f1.rename(name2);
    f2.rename(name1);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

RASModelVariables::RASModelVariables
(
    const fvMesh& mesh,
    const solverControl& SolverControl
)
:
    mesh_(mesh),
    solverControl_(SolverControl),

    TMVar1BaseName_(),
    TMVar2BaseName_(),
    nutBaseName_("nut"),

    TMVar1Ptr_(nullptr),
    TMVar2Ptr_(nullptr),
    nutPtr_(nullptr),
    distPtr_(nullptr),

    TMVar1InitPtr_(0),
    TMVar2InitPtr_(0),
    nutInitPtr_(0),

    TMVar1MeanPtr_(nullptr),
    TMVar2MeanPtr_(nullptr),
    nutMeanPtr_(nullptr)
{}


RASModelVariables::RASModelVariables
(
    const RASModelVariables& rmv
)
:
    mesh_(rmv.mesh_),
    solverControl_(rmv.solverControl_),

    TMVar1BaseName_(rmv.TMVar1BaseName_),
    TMVar2BaseName_(rmv.TMVar2BaseName_),
    nutBaseName_(rmv.nutBaseName_),
    TMVar1InitPtr_(rmv.TMVar1InitPtr_.size()),
    TMVar2InitPtr_(rmv.TMVar2InitPtr_.size()),
    nutInitPtr_(rmv.nutInitPtr_.size()),
    TMVar1MeanPtr_
    (
        variablesSet::allocateRenamedField(rmv.TMVar1MeanPtr_)
    ),
    TMVar2MeanPtr_
    (
        variablesSet::allocateRenamedField(rmv.TMVar2MeanPtr_)
    ),
    nutMeanPtr_
    (
        variablesSet::allocateRenamedField(rmv.nutMeanPtr_)
    )
{
    forAll(TMVar1InitPtr_, i)
    {
        if (rmv.TMVar1InitPtr_.get(i))
        {
            TMVar1InitPtr_.set
            (
                i,
                variablesSet::allocateRenamedField(rmv.TMVar1InitPtr_[i])
            );
        }
    }
    forAll(TMVar2InitPtr_, i)
    {
        if (rmv.TMVar2InitPtr_.get(i))
        {
            TMVar2InitPtr_.set
            (
                i,
                variablesSet::allocateRenamedField(rmv.TMVar2InitPtr_[i])
            );
        }
    }
    forAll(nutInitPtr_, i)
    {
        if (rmv.nutInitPtr_.get(i))
        {
            nutInitPtr_.set
            (
                i,
                variablesSet::allocateRenamedField(rmv.nutInitPtr_[i])
            );
        }
    }
}


autoPtr<RASModelVariables> RASModelVariables::clone() const
{
    return autoPtr<RASModelVariables>::New(*this);
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<RASModelVariables> RASModelVariables::New
(
    const fvMesh& mesh,
    const solverControl& SolverControl
)
{
    const IOdictionary modelDict
    (
        IOobject
        (
            turbulenceModel::propertiesName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false // Do not register
        )
    );

    word modelType("laminar"); // default to laminar

    const dictionary* dictptr = modelDict.findDict("RAS");

    if (dictptr)
    {
        // "RASModel" for v2006 and earlier
        dictptr->readCompat("model", {{"RASModel", -2006}}, modelType);
    }
    else
    {
        dictptr = &dictionary::null;
    }

    Info<< "Creating references for RASModel variables : " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            *dictptr,
            "RASModelVariables",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<RASModelVariables>(ctorPtr(mesh, SolverControl));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> RASModelVariables::nutJacobianVar1
(
    const singlePhaseTransportModel& laminarTransport
) const
{
    WarningInFunction
        << "jutJacobianVar1 not implemented for the current turbulence model."
        << "Returning zero field" << endl;

    return tmp<volScalarField>::New
    (
        IOobject
        (
            "nutJacobianVar1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    );
}


tmp<volScalarField> RASModelVariables::nutJacobianVar2
(
    const singlePhaseTransportModel& laminarTransport
) const
{
    WarningInFunction
        << "nutJacobianVar2 not implemented for the current turbulence model."
        << "Returning zero field" << endl;

    return tmp<volScalarField>::New
    (
        IOobject
        (
            "nutJacobianVar2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    );
}


void RASModelVariables::restoreInitValues()
{
    if (solverControl_.storeInitValues())
    {
        if (hasTMVar1())
        {
            variablesSet::restoreFieldInitialization
            (
                TMVar1InitPtr_,
                TMVar1Inst()
            );
        }
        if (hasTMVar2())
        {
            variablesSet::restoreFieldInitialization
            (
                TMVar2InitPtr_,
                TMVar2Inst()
            );
        }
        if (hasNut())
        {
            variablesSet::restoreFieldInitialization
            (
                nutInitPtr_,
                nutRefInst()
            );
        }
    }
}


void RASModelVariables::resetMeanFields()
{
    if (solverControl_.average())
    {
        Info<< "Resetting mean turbulent fields to zero" << endl;

        // Reset fields to zero
        if (hasTMVar1())
        {
            TMVar1MeanPtr_.ref() ==
                dimensionedScalar(TMVar1Inst().dimensions(), Zero);
        }
        if (hasTMVar2())
        {
            TMVar2MeanPtr_.ref() ==
                dimensionedScalar(TMVar2Inst().dimensions(), Zero);
        }
        if (hasNut())
        {
            nutMeanPtr_.ref() ==
                dimensionedScalar(nutRefInst().dimensions(), Zero);
        }
    }
}


void RASModelVariables::computeMeanFields()
{
    if (solverControl_.doAverageIter())
    {
        const label iAverageIter = solverControl_.averageIter();
        const scalar avIter(iAverageIter);
        const scalar oneOverItP1 = 1./(avIter + 1);
        const scalar mult = avIter*oneOverItP1;

        if (hasTMVar1())
        {
            TMVar1MeanPtr_.ref() ==
                (TMVar1MeanPtr_()*mult + TMVar1Inst()*oneOverItP1);
        }
        if (hasTMVar2())
        {
            TMVar2MeanPtr_.ref() ==
                (TMVar2MeanPtr_()*mult + TMVar2Inst()*oneOverItP1);
        }
        if (hasNut())
        {
            nutMeanPtr_.ref() ==
                (nutMeanPtr_()*mult + nutRefInst()*oneOverItP1);
        }
    }
}


void RASModelVariables::computeMeanUnsteadyFields()
{
    if (solverControl_.doAverageTime())
    {
        const scalar dt = mesh_.time().deltaTValue();
        const scalar elapsedTime
            = mesh_.time().value() - solverControl_.averageStartTime();
        scalar oneOverItP1 = dt/(elapsedTime + dt);
        scalar mult = elapsedTime/(elapsedTime + dt);
        if (hasTMVar1())
        {
            TMVar1MeanPtr_.ref() ==
                TMVar1MeanPtr_()*mult + TMVar1Inst()*oneOverItP1;
        }
        if (hasTMVar2())
        {
            TMVar2MeanPtr_.ref() ==
                TMVar2MeanPtr_()*mult + TMVar2Inst()*oneOverItP1;
        }
        if (hasNut())
        {
            nutMeanPtr_.ref() == nutMeanPtr_()*mult + nutRefInst()*oneOverItP1;
        }
    }
}


tmp<volSymmTensorField> RASModelVariables::devReff
(
    const singlePhaseTransportModel& laminarTransport,
    const volVectorField& U
) const
{
    return tmp<volSymmTensorField>::New
    (
        IOobject
        (
            "devRhoReff",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
      - (laminarTransport.nu() + nutRef())*dev(twoSymm(fvc::grad(U)))
    );
}


void RASModelVariables::correctBoundaryConditions
(
    const incompressible::turbulenceModel& turbulence
)
{
    if (hasTMVar1())
    {
        TMVar1Inst().correctBoundaryConditions();
        if (solverControl_.average())
        {
            TMVar1MeanPtr_.ref().correctBoundaryConditions();
        }
    }

    if (hasTMVar2())
    {
        TMVar2Inst().correctBoundaryConditions();
        if (solverControl_.average())
        {
            TMVar2MeanPtr_.ref().correctBoundaryConditions();
        }
    }

    if (hasNut())
    {
        nutRefInst().correctBoundaryConditions();
        if (solverControl_.average())
        {
            nutMeanPtr_.ref().correctBoundaryConditions();
        }
    }
}


void RASModelVariables::transfer(RASModelVariables& rmv)
{
    if (rmv.hasTMVar1() && hasTMVar1())
    {
        copyAndRename(TMVar1Inst(), rmv.TMVar1Inst());
    }

    if (rmv.hasTMVar2() && hasTMVar2())
    {
        copyAndRename(TMVar2Inst(), rmv.TMVar2Inst());
    }

    if (rmv.hasNut() && hasNut())
    {
        copyAndRename(nutRefInst(), rmv.nutRefInst());
    }

    if (rmv.hasDist() && hasDist())
    {
        copyAndRename(d(), rmv.d());
    }
}


void RASModelVariables::setWriteOption(IOobject::writeOption w)
{
    if (hasTMVar1())
    {
        TMVar1Inst().writeOpt() = w;
    }
    if (hasTMVar2())
    {
        TMVar2Inst().writeOpt() = w;
    }
    if (hasNut())
    {
        nutRefInst().writeOpt() = w;
    }
    if (solverControl_.doAverageTime())
    {
        if (hasTMVar1())
        {
            TMVar1MeanPtr_->writeOpt() = w;
        }
        if (hasTMVar2())
        {
            TMVar2MeanPtr_->writeOpt() = w;
        }
        if (hasNut())
        {
            nutMeanPtr_->writeOpt() = w;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
