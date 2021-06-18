/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "momentumError.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcLaplacian.H"
#include "turbulenceModel.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "topoSetSource.H"
#include "addToRunTimeSelectionTable.H"
#include "topoSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(momentumError, 0);
    addToRunTimeSelectionTable(functionObject, momentumError, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::functionObjects::momentumError::divDevRhoReff()
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    const auto& U = lookupObject<volVectorField>(UName_);
    tmp<volVectorField> tU(U);

    {
        auto* turb = findObject<cmpTurbModel>
        (
            turbulenceModel::propertiesName
        );

        if (turb)
        {
            tmp<volScalarField> trho = turb->rho();
            tmp<volScalarField> tnuEff = turb->nuEff();

            if (zoneSubSetPtr_)
            {
                const fvMeshSubset& subSetMesh= zoneSubSetPtr_->subSetMesh();

                tU.reset
                (
                    new volVectorField
                    (
                        subSetMesh.interpolate(U, false)
                    )

                );
                trho.reset
                (
                    new volScalarField
                    (
                        subSetMesh.interpolate(turb->rho(), false)
                    )
                );
                tnuEff.reset
                (
                    new volScalarField
                    (
                        subSetMesh.interpolate(turb->nuEff()(), false)
                    )
                );
            }
            return tmp<volVectorField>::New
            (
                "divDevRhoReff",
                - fvc::div
                (
                    (trho()*tnuEff())
                    *dev2(T(fvc::grad(tU()))),
                    "div(((rho*nuEff)*dev2(T(grad(U)))))"
                )
                - fvc::laplacian
                (
                    trho()*tnuEff(),
                    tU(),
                    "laplacian(nuEff,U)"
                )
            );
        }
    }

    {
        const auto* turb = findObject<icoTurbModel>
        (
            turbulenceModel::propertiesName
        );

        if (turb)
        {
            tmp<volScalarField> tnuEff = turb->nuEff();

            if (zoneSubSetPtr_)
            {
                const fvMeshSubset& subSetMesh= zoneSubSetPtr_->subSetMesh();

                tnuEff.reset
                (
                    new volScalarField
                    (
                        subSetMesh.interpolate(turb->nuEff()(), false)
                    )
                );
                tU.reset
                (
                    new volVectorField
                    (
                        subSetMesh.interpolate(U, false)
                    )

                );
            }

            return tmp<volVectorField>::New
            (
                "divDevRhoReff",
                - fvc::div
                (
                    tnuEff()*dev2(T(fvc::grad(tU()))),
                    "div(((nuEff)*dev2(T(grad(U)))))"
                )
                - fvc::laplacian(tnuEff(), tU(), "laplacian(nuEff,U)")
            );
        }
     }

     return volVectorField::null();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::momentumError::momentumError
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    pName_("p"),
    UName_("U"),
    phiName_("phi")
{
    read(dict);

    const auto& phi =lookupObject<surfaceScalarField>(phiName_);

    word momName(scopedName("momentError"));

    if (zoneSubSetPtr_)
    {
        const fvMesh& subMesh = zoneSubSetPtr_->subSetMesh().subMesh();

        auto* momentMapPtr
        (
            new volVectorField
            (
                IOobject
                (
                    scopedName("momentErrorMap"),
                    subMesh.time().timeName(),
                    subMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                subMesh,
                dimensionedVector(phi.dimensions()*dimVelocity/dimVolume, Zero)
            )
        );

        subMesh.objectRegistry::store(momentMapPtr);

        momName = scopedName("momentErrorZone");
    }

    auto* momentPtr
    (
        new volVectorField
        (
            IOobject
            (
                momName,
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector(phi.dimensions()*dimVelocity/dimVolume, Zero)
        )
    );

    mesh_.objectRegistry::store(momentPtr);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::momentumError::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {
        Info<< type() << " " << name() << ":" << nl;

        // Optional field name entries
        if (dict.readIfPresent<word>("p", pName_))
        {
            Info<< "    p: " << pName_ << endl;
        }
        if (dict.readIfPresent<word>("U", UName_))
        {
            Info<< "    U: " << UName_ << endl;
        }

        if (dict.readIfPresent<word>("phi", phiName_))
        {
            Info<< "    phi: " << phiName_ << endl;
        }

        if (dict.found("cellZones"))
        {
            zoneSubSetPtr_.reset(new zoneSubSet(name(), mesh_, dict));
        }

        return true;
    }

    return false;
}


void Foam::functionObjects::momentumError::calcMomentError()
{
    const auto& p = lookupObject<volScalarField>(pName_);
    const auto& U = lookupObject<volVectorField>(UName_);
    const auto& phi = lookupObject<surfaceScalarField>(phiName_);

    if (zoneSubSetPtr_)
    {
        zoneSubSetPtr_->subSetMesh().subMesh().fvSchemes::readOpt() =
            mesh_.fvSchemes::readOpt();

        zoneSubSetPtr_->subSetMesh().subMesh().fvSchemes::read();

        const fvMeshSubset& subFvMesh= zoneSubSetPtr_->subSetMesh();

        auto& momentErrMap =
            subFvMesh.subMesh().lookupObjectRef<volVectorField>
            (
                scopedName("momentErrorMap")
            );

        tmp<volScalarField> tp = subFvMesh.interpolate(p, false);
        tmp<volVectorField> tU = subFvMesh.interpolate(U, false);
        tmp<surfaceScalarField> tphi =
            subFvMesh.interpolate(phi, false);

         momentErrMap =
         (
            divDevRhoReff()
          + fvc::div(tphi, tU, "div(phi,U)")
          + fvc::grad(tp, "grad(p)")
         );
    }
    else
    {
        auto& momentErr =
            lookupObjectRef<volVectorField>(scopedName("momentError"));

        momentErr =  fvc::div(phi, U) + fvc::grad(p) + divDevRhoReff();
    }
}


bool Foam::functionObjects::momentumError::execute()
{
    calcMomentError();

    return true;
}


bool Foam::functionObjects::momentumError::write()
{
    Log << "    functionObjects::" << type() << " " << name();

    if (!zoneSubSetPtr_)
    {
        Log << " writing field: " << scopedName("momentError") << endl;

        const auto& momentErr =
            lookupObjectRef<volVectorField>(scopedName("momentError"));

        momentErr.write();
    }
    else
    {
        Log << " writing field: " << scopedName("momentErrorMap") << endl;

        const fvMeshSubset& subSetMesh= zoneSubSetPtr_->subSetMesh();
        const fvMesh& subMesh = subSetMesh.subMesh();

        const auto& momentErrMap =
            subMesh.lookupObject<volVectorField>
            (
                scopedName("momentErrorMap")
            );

        tmp<volVectorField> mapMomErr =
            zoneSubSetPtr_->mapToZone<vector>(momentErrMap);

        mapMomErr().write();
    }

    return true;
}


// ************************************************************************* //
