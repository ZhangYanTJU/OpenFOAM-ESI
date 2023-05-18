/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Louis Vittoz, SimScale GmbH
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

#include "fanMomentumSource.H"
#include "fvMatrices.H"
#include "DimensionedField.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"
#include "TableFile.H"
#include "turbulenceModel.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(fanMomentumSource, 0);
    addToRunTimeSelectionTable(option, fanMomentumSource, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions   * * * * * * * * * * * //

void Foam::fv::fanMomentumSource::writeProps
(
    const scalar gradP,
    const scalar flowRate
) const
{
    // Only write on output time
    if (mesh_.time().writeTime())
    {
        IOdictionary propsDict
        (
            IOobject
            (
                name_ + "Properties",
                mesh_.time().timeName(),
                "uniform",
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );
        propsDict.add("gradient", gradP);
        propsDict.add("flow_rate", flowRate);
        propsDict.regIOobject::write();
    }
}


void Foam::fv::fanMomentumSource::initializeUpstreamFaces()
{
    // First calculate the centre of gravity of the cell zone
    const vectorField& cellCentres = mesh_.cellCentres();
    const scalarField& cellVolumes = mesh_.cellVolumes();

    vector centreGravityCellZone = vector::zero;
    scalar cellZoneVolume = 0.;
    forAll(cells_, i)
    {
        const label cellI = cells_[i];
        const scalar cellVolume = cellVolumes[cellI];
        centreGravityCellZone += cellCentres[cellI]*cellVolume;
        cellZoneVolume += cellVolume;
    }

    reduce(centreGravityCellZone, sumOp<vector>());
    reduce(cellZoneVolume, sumOp<scalar>());

    centreGravityCellZone /= max(cellZoneVolume, SMALL);

    // Collect faces upstream of the centre of gavity
    const faceZone& fZone = mesh_.faceZones()[surroundingFaceZoneID_];
    const vectorField& faceCentreGravity = mesh_.faceCentres();

    upstreamFaceIDs_.setSize(fZone.size());
    upstreamFacePatchIDs_.setSize(fZone.size());

    label count = 0;
    forAll(fZone, i)
    {
        const label faceI = fZone[i];

        if
        (
            (flowDir_ & (faceCentreGravity[faceI] - centreGravityCellZone)) < 0.
        )
        {
            label faceId = -1;
            label facePatchId = -1;
            if (mesh_.isInternalFace(faceI))
            {
                faceId = faceI;
                facePatchId = -1;
            }
            else
            {
                facePatchId = mesh_.boundaryMesh().whichPatch(faceI);
                const polyPatch& pp = mesh_.boundaryMesh()[facePatchId];
                const auto* cpp = isA<coupledPolyPatch>(pp);

                if (cpp)
                {
                    faceId = cpp->owner() ? pp.whichFace(faceI) : -1;
                }
                else if (!isA<emptyPolyPatch>(pp))
                {
                    faceId = pp.whichFace(faceI);
                }
                else
                {
                    faceId = -1;
                    facePatchId = -1;
                }
            }
            if (faceId >= 0)
            {
                upstreamFacePatchIDs_[count] = facePatchId;
                upstreamFaceIDs_[count] = faceId;
                count++;
            }
        }
    }

    upstreamFaceIDs_.setSize(count);
    upstreamFacePatchIDs_.setSize(count);

    // Fill cellsInZones_ with all cell IDs
    forAll(cells_, i)
    {
        cellsInZones_.insert(cells_[i]);
    }

    // Sanity check
    const labelUList& owners = mesh_.owner();
    const labelUList& neighbours = mesh_.neighbour();
    forAll(upstreamFaceIDs_, i)
    {
        if (upstreamFacePatchIDs_[i] == -1)
        {
            const label faceI = upstreamFaceIDs_[i];
            const label own = owners[faceI];
            const label nei = neighbours[faceI];
            
            // To be valid: one cell has to be inside the cellZone and the other
            // one, outside
            if (cellsInZones_.found(own) == cellsInZones_.found(nei))
            {
                FatalErrorInFunction
                    << "It seems that the faceZone is not part of the cellZone "
                    << "boundaries."
                    << abort(FatalError);
            }
        }
    }
}


template<typename flowRateFunctorPatch, typename flowRateFunctor>
Foam::scalar Foam::fv::fanMomentumSource::calculateFlowRate
(
    flowRateFunctorPatch fPatch,
    flowRateFunctor f
) const
{
    // Calculate the flow rate over the upstream faces
    scalarList phif(upstreamFaceIDs_.size());

    const labelUList& owners = mesh_.owner();

    forAll(upstreamFaceIDs_, i)
    {
        const label faceI = upstreamFaceIDs_[i];
        if (upstreamFacePatchIDs_[i] != -1)
        {
            const label patchI = upstreamFacePatchIDs_[i];
            phif[i] = fPatch(patchI,faceI);
        }
        else
        {
            const label own = owners[faceI];
            if (cellsInZones_.found(own))
            {
                // Owner is in cell zone, which means that neighbour is not.
                // Positive flux means that the flux is going out of the domain.
                // Flip the sign.
                phif[i] = -f(faceI);
            }
            else // if (cellsInZones_.found(nei))
            {
                // Neighbour is in cell zone, which means that owner is not.
                // Positive flux means that the flux is coming into the domain.
                // Don't flip the sign.
                phif[i] = f(faceI);
            }
        }
    }

    return gSum(phif);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::fanMomentumSource::fanMomentumSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(sourceName, modelType, dict, mesh),
    fanCurve_(Function1<scalar>::New("fanCurve", dict)),
    flowDir_(coeffs_.get<vector>("flowDir")),
    thickness_(coeffs_.get<scalar>("thickness")),
    gradPFan_(0.0),
    upstreamFaceIDs_(),
    upstreamFacePatchIDs_(),
    rho_(nullptr)
{
    // Skip all the checks if the source term has been deactivated
    // because there are no selected cells.
    // Such a situation typically occurs for multiple fluid regions
    if (fv::option::isActive())
    {
        const word faceZoneName = coeffs_.get<word>("faceZone");

        surroundingFaceZoneID_ = mesh_.faceZones().findZoneID(faceZoneName);

        if (surroundingFaceZoneID_ < 0)
        {
            FatalErrorInFunction
                << type() << " " << this->name() << ": "
                << "    Unknown face zone name: " << faceZoneName
                << ". Valid face zones are: " << mesh_.faceZones().names()
                << exit(FatalError);
        }

        if (mag(flowDir_) < SMALL)
        {
            FatalErrorInFunction
                << "Detected zero-vector for flowDir. Check your settings."
                << exit(FatalError);
        }
        else
        {
            flowDir_ /= mag(flowDir_);
        }

        if (thickness_ <= 0.)
        {
            FatalErrorInFunction
                << "The thickness of the fan model region cannot be negative"
                << " or null."
                << exit(FatalError);
        }


        coeffs_.readEntry("fields", fieldNames_);

        if (fieldNames_.size() != 1)
        {
            FatalErrorInFunction
                << "Source can only be applied to a single field. Current "
                << "settings are:" << fieldNames_ << exit(FatalError);
        }

        fv::option::resetApplied();

        // Read the initial pressure gradient from file if it exists
        IFstream propsFile
        (
            mesh_.time().timePath()/"uniform"/(name_ + "Properties")
        );

        if (propsFile.good())
        {
            Info<< "    Reading pressure gradient from file" << endl;
            dictionary propsDict(propsFile);
            propsDict.readEntry("gradient", gradPFan_);
        }

        Info<< "    Initial pressure gradient = " << gradPFan_ << nl << endl;

        if (coeffs_.found("rho"))
        {
            rho_.reset(new scalar(coeffs_.get<scalar>("rho")));
        }

        initializeUpstreamFaces();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::fanMomentumSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    volVectorField::Internal Su
    (
        IOobject
        (
            name_ + fieldNames_[fieldi] + "Sup",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(eqn.dimensions()/dimVolume, Zero)
    );

    const auto& phi = mesh().lookupObject<surfaceScalarField>("phi");

    if (phi.dimensions() != dimVelocity*dimArea)
    {
        FatalErrorInFunction
            << "You called incompressible variant of addSup for case with "
            << "a mass flux and not volumetric flux. This is not allowed."
            << abort(FatalError);
    }

    if (!rho_.good() || rho_.ref() < VSMALL)
    {
        FatalErrorInFunction
            << "You called incompressible addSup without or with "
            << "zero value reference density."
            << abort(FatalError);
    }

    // Lambda functions passed as arguments to calculateFlowRate
    auto getVolumetricFlowRatePatch =
        [&phi](const label patchI, const label faceI)
        {
            return phi.boundaryField()[patchI][faceI];
        };
    auto getVolumetricFlowRate =
        [&phi](const label faceI)
        {
            return phi[faceI];
        };

    const scalar flowRate =
        calculateFlowRate(getVolumetricFlowRatePatch, getVolumetricFlowRate);

    // Pressure drop for this flow rate
    // if flow rate is negative, pressure is clipped at the static pressure
    gradPFan_ = fanCurve_->value(max(flowRate, scalar(0)))/thickness_/rho_.ref();

    // Create the source term
    UIndirectList<vector>(Su, cells_) = flowDir_*gradPFan_;

    eqn += Su;

    writeProps(gradPFan_, flowRate);
}


void Foam::fv::fanMomentumSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    volVectorField::Internal Su
    (
        IOobject
        (
            name_ + fieldNames_[fieldi] + "Sup",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(eqn.dimensions()/dimVolume, Zero)
    );

    const auto& phi = mesh().lookupObject<surfaceScalarField>("phi");

    if (phi.dimensions() != dimVelocity*dimArea*dimDensity)
    {
        FatalErrorInFunction
            << "You called compressible variant of addSup for case with "
            << "a volumetric flux and not mass flux. This is not allowed."
            << abort(FatalError);
    }

    const surfaceScalarField rhof = fvc::interpolate(rho);

    // Lambda functions passed as arguments to calculateFlowRate
    auto getVolumetricFlowRatePatch =
        [&phi, &rhof](const label patchI, const label faceI)
        {
            return phi.boundaryField()[patchI][faceI]
                  /rhof.boundaryField()[patchI][faceI];
        };
    auto getVolumetricFlowRate =
        [&phi, &rhof](const label faceI)
        {
            return phi[faceI]/rhof.internalField()[faceI];
        };

    const scalar flowRate =
        calculateFlowRate(getVolumetricFlowRatePatch, getVolumetricFlowRate);

    // Pressure drop for this flow rate
    // if flow rate is negative, pressure is clipped at the static pressure
    gradPFan_ = fanCurve_->value(max(flowRate, scalar(0)))/thickness_;

    // Create the source term
    UIndirectList<vector>(Su, cells_) = flowDir_*gradPFan_;

    eqn += Su;

    writeProps(gradPFan_, flowRate);

}


bool Foam::fv::fanMomentumSource::read(const dictionary& dict)
{
    NotImplemented;

    return false;
}


// ************************************************************************* //
