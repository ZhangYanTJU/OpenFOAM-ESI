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
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
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
    for (const label celli : cells_)
    {
        const scalar cellVolume = cellVolumes[celli];
        centreGravityCellZone += cellCentres[celli]*cellVolume;
        cellZoneVolume += cellVolume;
    }

    reduce(centreGravityCellZone, sumOp<vector>());
    reduce(cellZoneVolume, sumOp<scalar>());

    centreGravityCellZone /= max(cellZoneVolume, SMALL);

    // Collect faces upstream of the centre of gavity
    const faceZone& fZone = mesh_.faceZones()[surroundingFaceZoneID_];
    const vectorField& faceCentreGravity = mesh_.faceCentres();

    upstreamPatchFaceInfo_.resize_nocopy(fZone.size());

    label count = 0;
    for (const label facei : fZone)
    {
        if
        (
            (flowDir_ & (faceCentreGravity[facei] - centreGravityCellZone)) < 0.
        )
        {
            labelPair patchFaceInfo(-1, -1);

            if (mesh_.isInternalFace(facei))
            {
                // Patch ID already set to -1, set only the face ID
                patchFaceInfo.second() = facei;
            }
            else
            {
                patchFaceInfo.first() = mesh_.boundaryMesh().whichPatch(facei);
                const polyPatch& pp = mesh_.boundaryMesh()[patchFaceInfo.first()];
                const auto* cpp = isA<coupledPolyPatch>(pp);

                if (cpp)
                {
                    patchFaceInfo.second() =
                        cpp->owner()
                      ? pp.whichFace(facei)
                      : -1;
                }
                else if (!isA<emptyPolyPatch>(pp))
                {
                    patchFaceInfo.second() = pp.whichFace(facei);
                }
                // else both face ID and patch ID remain at -1
            }

            // If this is an upstream face, set it in the list
            if (patchFaceInfo.second() >= 0)
            {
                upstreamPatchFaceInfo_[count] = patchFaceInfo;
                count++;
            }
        }
    }

    upstreamPatchFaceInfo_.setSize(count);

    // Fill cellsInZones_ with all cell IDs
    for (const label celli : cells_)
    {
        cellsInZones_.insert(celli);
    }

    // Sanity check
    const labelUList& owners = mesh_.owner();
    const labelUList& neighbours = mesh_.neighbour();
    for (const labelPair& patchFaceInfo : upstreamPatchFaceInfo_)
    {
        if (patchFaceInfo.first() == -1)
        {
            const label facei = patchFaceInfo.second();
            const label own = owners[facei];
            const label nei = neighbours[facei];

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


template<typename FlowRateFunctor>
Foam::scalar
Foam::fv::fanMomentumSource::calculateFlowRate(FlowRateFunctor f) const
{
    // Calculate the flow rate through the upstream faces
    scalarList phif(upstreamPatchFaceInfo_.size());

    const labelUList& owners = mesh_.owner();

    forAll(upstreamPatchFaceInfo_, i)
    {
        const labelPair& patchFaceInfo = upstreamPatchFaceInfo_[i];

        // Sign of the flux needs to be flipped if this is an internal face
        // whose owner is found in the cell zone
        phif[i] =
            patchFaceInfo.first() < 0
        && cellsInZones_.found(owners[patchFaceInfo.second()])
         ? -f(patchFaceInfo)
         : f(patchFaceInfo);
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
    upstreamPatchFaceInfo_(),
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

    // Lambda function passed as argument to calculateFlowRate
    auto getVolumetricFlowRate =
        [&phi](const labelPair& patchFaceInfo)
        {
            return
            (
                patchFaceInfo.first() < 0
              ? phi[patchFaceInfo.second()]
              : phi.boundaryField()[patchFaceInfo]
            );
        };

    const scalar flowRate = calculateFlowRate(getVolumetricFlowRate);

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

    if (phi.dimensions() != dimMass/dimTime)
    {
        FatalErrorInFunction
            << "You called compressible variant of addSup for case with "
            << "a volumetric flux and not mass flux. This is not allowed."
            << abort(FatalError);
    }

    const surfaceScalarField rhof = fvc::interpolate(rho);

    // Lambda function passed as argument to calculateFlowRate
    auto getVolumetricFlowRate =
        [&phi, &rhof](const labelPair& patchFaceInfo)
        {
            return
            (
                patchFaceInfo.first() < 0
              ? phi[patchFaceInfo.second()]/
                rhof.internalField()[patchFaceInfo.second()]
              : phi.boundaryField()[patchFaceInfo]/
                rhof.boundaryField()[patchFaceInfo]
            );
        };

    const scalar flowRate = calculateFlowRate(getVolumetricFlowRate);

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
