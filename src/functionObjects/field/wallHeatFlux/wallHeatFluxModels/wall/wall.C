/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenCFD Ltd.
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

#include "wall.H"
#include "turbulentFluidThermoModel.H"
#include "solidThermo.H"
#include "surfaceInterpolate.H"
#include "fvcSnGrad.H"
#include "wallPolyPatch.H"
#include "turbulentFluidThermoModel.H"
#include "multiphaseInterSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallHeatFluxModels
{
    defineTypeNameAndDebug(wall, 0);
    addToRunTimeSelectionTable
    (
        wallHeatFluxModel,
        wall,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::wallHeatFluxModels::wall::writeFileHeader(Ostream& os)
{
    writeHeader(os, "Wall heat-flux");
    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    writeTabbed(os, "integral");
    os  << endl;

    writtenHeader_ = true;
}


void Foam::wallHeatFluxModels::wall::calcHeatFlux
(
    const volScalarField& alpha,
    const volScalarField& he,
    volScalarField& wallHeatFlux
)
{
    volScalarField::Boundary& wallHeatFluxBf = wallHeatFlux.boundaryFieldRef();

    const volScalarField::Boundary& heBf = he.boundaryField();

    const volScalarField::Boundary& alphaBf = alpha.boundaryField();

    const labelHashSet& patches = patchSet();
    for (const label patchi : patches)
    {
        wallHeatFluxBf[patchi] = alphaBf[patchi]*heBf[patchi].snGrad();
    }


    const auto* qrPtr = mesh().cfindObject<volScalarField>(qrName());

    if (qrPtr)
    {
        const volScalarField::Boundary& radHeatFluxBf = qrPtr->boundaryField();

        for (const label patchi : patches)
        {
            wallHeatFluxBf[patchi] -= radHeatFluxBf[patchi];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallHeatFluxModels::wall::wall
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& name,
    const word objName,
    functionObjects::stateFunctionObject& state
)
:
    wallHeatFluxModel(dict, mesh, name, objName, state)
{
    auto* wallHeatFluxPtr
    (
        new volScalarField
        (
            IOobject
            (
                objName,
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar(dimMass/pow3(dimTime), Zero)
        )
    );

    mesh.objectRegistry::store(wallHeatFluxPtr);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::wallHeatFluxModels::wall::read(const dictionary& dict)
{
    if (!wallHeatFluxModel::read(dict))
    {
        return false;
    }


    qrName_ = dict.getOrDefault<word>("qr", "qr");


    Info<< state().type() << " " << state().name() << ":" << nl;

    patchSet_ = mesh().boundaryMesh().patchSet
    (
        dict.getOrDefault<wordRes>("patches", wordRes())
    );

    const polyBoundaryMesh& pbm = mesh().boundaryMesh();
    if (patchSet_.empty())
    {
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                patchSet_.insert(patchi);
            }
        }

        Info<< "    processing all wall patches" << nl << endl;
    }
    else
    {
        Info<< "    processing wall patches: " << nl;
        labelHashSet filteredPatchSet;
        for (const label patchi : patchSet_)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                filteredPatchSet.insert(patchi);
                Info<< "        " << pbm[patchi].name() << endl;
            }
            else
            {
                WarningInFunction
                    << "Requested wall heat-flux on non-wall boundary "
                    << "type patch: " << pbm[patchi].name() << endl;
            }
        }

        Info<< endl;

        patchSet_ = filteredPatchSet;
    }


    if (writeFile::canResetFile())
    {
        writeFile::resetFile(objName());
    }

    if (writeFile::canWriteHeader())
    {
        writeFileHeader(file());
    }


    return true;
}


bool Foam::wallHeatFluxModels::wall::execute()
{
    auto& wallHeatFlux = mesh().lookupObjectRef<volScalarField>(objName());

    if
    (
        mesh().foundObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        )
    )
    {
        const compressible::turbulenceModel& turbModel =
            mesh().lookupObject<compressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
            );

        calcHeatFlux
        (
            turbModel.alphaEff()(),
            turbModel.transport().he(),
            wallHeatFlux
        );
    }
    else if (mesh().foundObject<fluidThermo>(fluidThermo::dictName))
    {
        const fluidThermo& thermo =
            mesh().lookupObject<fluidThermo>(fluidThermo::dictName);

        calcHeatFlux
        (
            thermo.alpha(),
            thermo.he(),
            wallHeatFlux
        );
    }
    else if (mesh().foundObject<solidThermo>(solidThermo::dictName))
    {
        const solidThermo& thermo =
            mesh().lookupObject<solidThermo>(solidThermo::dictName);

        calcHeatFlux(thermo.alpha(), thermo.he(), wallHeatFlux);
    }
    else if
    (
        mesh().foundObject<multiphaseInterSystem>
            (multiphaseInterSystem::phasePropertiesName)
    )
    {
        const auto& thermo = mesh().lookupObject<multiphaseInterSystem>
        (
            multiphaseInterSystem::phasePropertiesName
        );

        calcHeatFlux(thermo.kappaEff()(), thermo.T(), wallHeatFlux);
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find compressible turbulence model in the "
            << "database" << exit(FatalError);
    }

    const fvPatchList& patches = mesh().boundary();

    const surfaceScalarField::Boundary& magSf = mesh().magSf().boundaryField();

    const labelHashSet& patchset = patchSet();
    for (const label patchi : patchset)
    {
        const fvPatch& pp = patches[patchi];

        const scalarField& hfp = wallHeatFlux.boundaryField()[patchi];

        const scalar minHfp = gMin(hfp);
        const scalar maxHfp = gMax(hfp);
        const scalar integralHfp = gSum(magSf[patchi]*hfp);

        if (Pstream::master())
        {
            writeCurrentTime(file());

            file()
                << token::TAB << pp.name()
                << token::TAB << minHfp
                << token::TAB << maxHfp
                << token::TAB << integralHfp
                << endl;
        }

        if (state().log)
        {
            Info<< "    min/max/integ(" << pp.name() << ") = "
                << minHfp << ", " << maxHfp << ", " << integralHfp << endl;
        }

        state().setResult("min(" + pp.name() + ")", minHfp);
        state().setResult("max(" + pp.name() + ")", maxHfp);
        state().setResult("int(" + pp.name() + ")", integralHfp);
    }


    return true;
}


bool Foam::wallHeatFluxModels::wall::write()
{
    const auto& wallHeatFlux =
        mesh().lookupObject<volScalarField>(objName());

    if (state().log)
    {
        Info<< "    writing field " << wallHeatFlux.name() << endl;
    }

    wallHeatFlux.write();

    return true;
}


// ************************************************************************* //
