/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2022 PCOpt/NTUA
    Copyright (C) 2013-2022 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "incompressibleAdjointVars.H"
#include "adjointBoundaryCondition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(incompressibleAdjointVars, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

incompressibleAdjointVars::incompressibleAdjointVars
(
    fvMesh& mesh,
    solverControl& SolverControl,
    objectiveManager& objManager,
    incompressibleVars& primalVars
)
:
    incompressibleAdjointMeanFlowVars(mesh, SolverControl, primalVars),
    objectiveManager_(objManager),

    adjointTurbulence_
    (
        incompressibleAdjoint::adjointRASModel::New
        (
            primalVars_,
            *this,
            objManager
        )
    ),
    writeFields_
    (
        SolverControl.solverDict().getOrDefault<bool>("writeFields", true)
    )
{
    if (!writeFields_)
    {
        setWriteOption(IOobject::NO_WRITE);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void incompressibleAdjointVars::restoreInitValues()
{
    if (solverControl_.storeInitValues())
    {
        Info<< "Restoring field values to initial ones" << endl;
        variablesSet::resetAdjointField(paInst());
        variablesSet::resetAdjointField(UaInst());
        variablesSet::resetAdjointField(phiaInst());
        adjointTurbulence_().restoreInitValues();
    }
}


void incompressibleAdjointVars::resetMeanFields()
{
    if (solverControl_.average())
    {
        Info<< "Resetting adjoint mean fields to zero" << endl;

        // Reset fields to zero
        paMeanPtr_() == dimensionedScalar(paPtr_().dimensions(), Zero);
        UaMeanPtr_() == dimensionedVector(UaPtr_().dimensions(), Zero);
        phiaMeanPtr_() == dimensionedScalar(phiaPtr_().dimensions(), Zero);
        adjointTurbulence_().resetMeanFields();

        // Reset averaging iteration index to 0
        solverControl_.averageIter() = 0;
    }
}


void incompressibleAdjointVars::computeMeanFields()
{
    if (solverControl_.doAverageIter())
    {
        Info<< "Averaging adjoint fields" << endl;
        label& iAverageIter = solverControl_.averageIter();
        scalar avIter(iAverageIter);
        scalar oneOverItP1 = 1./(avIter+1);
        scalar mult = avIter*oneOverItP1;
        paMeanPtr_() == paMeanPtr_()  *mult + paPtr_()  *oneOverItP1;
        UaMeanPtr_() == UaMeanPtr_()  *mult + UaPtr_()  *oneOverItP1;
        phiaMeanPtr_() == phiaMeanPtr_()*mult + phiaPtr_()*oneOverItP1;
        adjointTurbulence_().computeMeanFields();
        ++iAverageIter;
    }
}


void incompressibleAdjointVars::computeMeanUnsteadyFields
(
    const scalar endTime
)
{
    if (solverControl_.doAverageTime())
    {
        Info<< "Averaging fields" << endl;
        const scalar dt = mesh_.time().deltaTValue();
        const scalar elapsedTime = endTime - mesh_.time().value();
        scalar oneOverItP1 = dt/(elapsedTime + dt);
        scalar mult = elapsedTime/(elapsedTime + dt);
        paMeanPtr_() == paMeanPtr_()*mult + paInst()*oneOverItP1;
        UaMeanPtr_() == UaMeanPtr_()*mult + UaInst()*oneOverItP1;
        phiaMeanPtr_() == phiaMeanPtr_()*mult + phiaInst()*oneOverItP1;
        adjointTurbulence_().computeMeanUnsteadyFields(endTime);
    }
}


void incompressibleAdjointVars::nullify()
{
    incompressibleAdjointMeanFlowVars::nullify();
    adjointTurbulence_->nullify();
}


void incompressibleAdjointVars::updatePrimalBasedQuantities()
{
    /*
    // WIP
    for (fvPatchVectorField& pf : UaInst().boundaryFieldRef())
    {
        if (isA<adjointBoundaryCondition<vector>>(pf))
        {
            adjointBoundaryCondition<vector>& adjointBC =
                refCast<adjointBoundaryCondition<vector>>(pf);
            adjointBC.updatePrimalBasedQuantities();
        }
    }

    for (fvPatchScalarField& pf : paInst().boundaryFieldRef())
    {
        if (isA<adjointBoundaryCondition<scalar>>(pf))
        {
            adjointBoundaryCondition<scalar>& adjointBC =
                refCast<adjointBoundaryCondition<scalar>>(pf);
            adjointBC.updatePrimalBasedQuantities();
        }
    }
    */
}


void incompressibleAdjointVars::write()
{
    if (writeFields_)
    {
        incompressibleAdjointMeanFlowVars::write();
        adjointTurbulence_->writeFields();
    }
}


void incompressibleAdjointVars::setWriteOption(IOobject::writeOption w)
{
    if (!writeFields_)
    {
        w = IOobject::NO_WRITE;
    }

    incompressibleAdjointMeanFlowVars::setWriteOption(w);
    adjointTurbulence_->setWriteOption(w);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
