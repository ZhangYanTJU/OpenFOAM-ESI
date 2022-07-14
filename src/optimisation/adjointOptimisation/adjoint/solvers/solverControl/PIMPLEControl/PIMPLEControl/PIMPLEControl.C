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

#include "PIMPLEControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PIMPLEControl, 0);
    defineRunTimeSelectionTable(PIMPLEControl, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PIMPLEControl::PIMPLEControl
(
    fvMesh& mesh,
    const word& managerType,
    const solver& solver
)
:
    solverControl(solver),
    pimpleControl(mesh, "PIMPLE", false),
    managerType_(managerType),
    pRefCell_(Zero),
    pRefValue_(Zero),
    timeManip_
    (
        const_cast<unsteadyTimeManipulation&>
        (
            unsteadyTimeManipulation::New(mesh)
        )
    )
{
    read();
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::PIMPLEControl> Foam::PIMPLEControl::New
(
    fvMesh& mesh,
    const word& managerType,
    const solver& solver
)
{
    auto* ctorPtr = dictionaryConstructorTable(managerType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "control",
            managerType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<PIMPLEControl>(ctorPtr(mesh, managerType, solver));
}


bool Foam::PIMPLEControl::read()
{
    pimpleControl::read();
    solverControl::read();

    // Always store initial field values in unsteady runs
    storeInitValues_ = true;

    /*Info<< nl << algorithmName_;

    if (nCorrPIMPLE_ > 1)
    {
        if (residualControl_.empty())
        {
            Info<< ": no residual control data found. "
                << "Calculations will employ " << nCorrPIMPLE_
                << " corrector loops" << nl;
        }
        else
        {
            Info<< ": max iterations = " << nCorrPIMPLE_ << nl;

            for (const fieldData& ctrl : residualControl_)
            {
                Info<< "    field " << ctrl.name << token::TAB
                    << ": relTol " << ctrl.relTol
                    << ", tolerance " << ctrl.absTol
                    << nl;
            }
        }
    }
    else
    {
        Info<< ": Operating solver in PISO mode" << nl;
    }

    Info<< endl;*/

    return true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PIMPLEControl::updateAverageStartTime()
{
    averageStartTime_ += timeManip_.span();
    DebugInfo
        << "Set averageStartTime to " << averageStartTime_ << endl;
}


// ************************************************************************* //
