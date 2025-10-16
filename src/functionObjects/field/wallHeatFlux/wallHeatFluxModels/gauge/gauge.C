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

#include "gauge.H"
#include "physicoChemicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallHeatFluxModels
{
    defineTypeNameAndDebug(gauge, 0);
    addToRunTimeSelectionTable
    (
        wallHeatFluxModel,
        gauge,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::wallHeatFluxModels::gauge::calcRadiantExitance() const
{
    // Radiant exitance (Wikipedia EN): https://w.wiki/CdZz
    return e_*constant::physicoChemical::sigma.value()*Foam::pow4(Tgauge_);
}


Foam::tmp<Foam::labelField>
Foam::wallHeatFluxModels::gauge::identifyProbeCells() const
{
    // Fetch the neighbour cells of the operand patch
    const polyPatch& patch = mesh().boundaryMesh()[patchID_];
    const labelUList& faceCells = patch.faceCells();
    const label patchFaceStart = patch.start();

    // Fetch the patch faces that correspond to the specified probes
    const labelList& probeFaces = patchProber::faces();

    // Initialise the size of the operand cells
    auto tprobeCells = tmp<labelField>::New(probeFaces.size(), -1);
    labelField& probeCells = tprobeCells.ref();

    // Store the indices of the cells that contain the probed patch faces
    forAll(probeFaces, facei)
    {
        if (activeFaces_[facei])
        {
            const label patchFaceGlobali = probeFaces[facei];
            const label patchFaceLocali = patchFaceGlobali - patchFaceStart;

            probeCells[facei] = faceCells[patchFaceLocali];
        }
    }

    return tprobeCells;
}


bool Foam::wallHeatFluxModels::gauge::writeFileHeader(Ostream& os)
{
    writeHeader(os, typeName);
    os  << "# Patch," << mesh().boundaryMesh()[patchID_].name() << nl
        << "# Tgauge," << Tgauge_ << nl
        << "# Absorptivity," << a_ << nl
        << "# Emissivity," << e_ << nl
        << "# ProbeID,Original Location,Patch Location" << nl;

    const pointField& oldPoints = patchProber::oldPoints();
    const pointField& probeLocs = patchProber::probeLocations();

    for (label i = 0; i < szProbes_; ++i)
    {
        const vector& oldPoint = oldPoints[i];
        const vector& probeLoc = probeLocs[i];

        os  << '#' << ' ' << i
            << ',' << oldPoint.x() << ',' << oldPoint.y() << ',' << oldPoint.z()
            << ',' << probeLoc.x() << ',' << probeLoc.y() << ',' << probeLoc.z()
            << nl;
    }

    os  << "# Time";
    for (int i = 0; i < szProbes_; ++i)
    {
        os  << ',' << i;
    }

    os  << endl;
    return true;
}


bool Foam::wallHeatFluxModels::gauge::calcConvectiveHeatFlux()
{
    // Fetch the turbulent thermal diffusivity field
    const auto* alphatPtr = mesh().cfindObject<volScalarField>(alphatName_);

#ifdef FULLDEBUG
    if (!alphatPtr)
    {
        FatalErrorInFunction
            << "Turbulent thermal diffusivity field, alphat, is not available."
            << nl << "alphat: " << alphatName_
            << exit(FatalError);
    }
#endif

    // Compute the term: C_p*(alpha + alpha_t)
    tmp<volScalarField> tkappaEff = thermo_.kappaEff(*alphatPtr);
    const volScalarField& kappaEff = tkappaEff.cref();

    // Sample the term above at patch-probe locations
    tmp<scalarField> tsampledKappaEff = patchProber::sample(kappaEff);
    const scalarField& sampledKappaEff = tsampledKappaEff.cref();


    // Compute the patch-normal gradient of temperature with respect to the
    // specified gauge temperature
    tmp<scalarField> tdTdn = calcdTdn();
    const scalarField& dTdn = tdTdn.cref();

#ifdef FULLDEBUG
    if (sampledKappaEff.size() != dTdn.size())
    {
        FatalErrorInFunction
            << "Size mismatch: kappaEff samples (" << sampledKappaEff.size()
            << ") vs dTdn (" << dTdn.size() << ")." << nl
            << "Probe selection/order must match." << nl
            << exit(FatalError);
    }
#endif

    // Compute q''_conv = k_eff * dT/dn
    tmp<scalarField> tq = sampledKappaEff*dTdn;
    const scalarField& q = tq.cref();

    // Be picky - Ensure storage exists and matches size
    if (qConvPtr_ && (qConvPtr_->size() == q.size()))
    {
        *qConvPtr_ = q;
    }
    else
    {
        IOobject io
        (
            "qConv",
            mesh().time().timeName(),
            mesh()
        );

        qConvPtr_.reset(new scalarIOField(io, q));
    }


    // Note that for inactive faces q=0 since dT/dn=0
    Pstream::listReduce(*qConvPtr_, sumOp<scalar>());

    return true;
}


Foam::tmp<Foam::scalarField>
Foam::wallHeatFluxModels::gauge::calcdTdn() const
{
    // Allocate and initialise the storage for dT/dn
    auto tdTdn = tmp<scalarField>::New(szProbes_, 0);
    scalarField& dTdn = tdTdn.ref();

    // Fetch the temperature field
    const auto* Tptr = mesh().cfindObject<volScalarField>(Tname_);

#ifdef FULLDEBUG
    if (!Tptr)
    {
        FatalErrorInFunction
            << "Temperature field, T, is not available."
            << nl << "T: " << Tname_
            << exit(FatalError);
    }
#endif

    // Fetch the internal fields for temperature and distance-to-patch fields
    const scalarField& Ti = Tptr->primitiveField();
    const auto& pdc = mesh().deltaCoeffs().boundaryField()[patchID_];


    // Compute dT/dn for active faces and their owner cells
    forAll(cells_, probei)
    {
        if (activeFaces_[probei])
        {
            const label cellID = cells_[probei];
            dTdn[probei] = (Ti[cellID] - Tgauge_)/max(pdc[probei], SMALL);
        }
    }

    return tdTdn;
}


bool Foam::wallHeatFluxModels::gauge::calcRadiativeHeatFlux()
{
    // Fetch the incident radiative heat-flux field
    const auto* qinPtr = mesh().getObjectPtr<volScalarField>(qinName_);

#ifdef FULLDEBUG
    if (!qinPtr)
    {
        FatalErrorInFunction
            << "Incident radiative heat-flux field, qin, is not available."
            << nl << "qin: " << qinName_
            << exit(FatalError);
    }
#endif

    const volScalarField& qin = *qinPtr;

    // Sample the incident radiative heat-flux field at patch-probe locations
    tmp<scalarField> tq  = patchProber::sample(qin);
    const scalarField& q = tq.cref();

    // Compute q''_rad = a * q_in - M_e
    tmp<scalarField> tqrad = a_*q - Me_;
    scalarField& qrad = tqrad.ref();

    // Overwrite q''_rad=0 for inactive probed faces
    forAll(qrad, probei)
    {
        if (!activeFaces_[probei])
        {
            qrad[probei] = 0;
        }
    }


    // Be picky - Ensure storage exists and matches size
    if (qRadPtr_ && (qRadPtr_->size() == qrad.size()))
    {
        *qRadPtr_ = qrad;
    }
    else
    {
        IOobject io
        (
            "qRad",
            mesh().time().timeName(),
            mesh()
        );

        qRadPtr_.reset(new scalarIOField(io, qrad));
    }

    // Note that for inactive faces qrad=0
    Pstream::listReduce(*qRadPtr_, sumOp<scalar>());

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallHeatFluxModels::gauge::gauge
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& name,
    const word objName,
    functionObjects::stateFunctionObject& state
)
:
    wallHeatFluxModel(dict, mesh, name, objName, state),
    patchProber(mesh, dict),
    thermo_(mesh.lookupObject<fluidThermo>(fluidThermo::dictName))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::wallHeatFluxModels::gauge::read(const dictionary& dict)
{
    if (!wallHeatFluxModel::read(dict) || !patchProber::read(dict))
    {
        return false;
    }


    // Do not proceed if the operand heat-flux field is not registered
    qinName_ = dict.getOrDefault<word>("qin", "qin");

    const auto* qinPtr = mesh().cfindObject<volScalarField>(qinName_);
    if (!qinPtr)
    {
        FatalErrorInFunction
            << "Incident radiative heat-flux field, qin, is not available."
            << nl << "qin: " << qinName_
            << exit(FatalError);
    }


    // Fetch the patch-prober data
    const labelList& patchIDs = patchProber::patchIDs();
    if (patchIDs.size() == 1)
    {
        patchID_ = patchIDs[0];
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "This function object can operate on only a single patch." << nl
            << "The number of input patches: " << patchIDs.size()
            << exit(FatalIOError);
    }


    // Fetch the number of probes, and verify the consistent sizing
    szProbes_ = patchProber::probeLocations().size();

    const labelList& probeFaces = patchProber::faces();

    if (szProbes_ != probeFaces.size())
    {
        FatalErrorInFunction
            << "Size mismatch: szProbes (" << szProbes_
            << ") vs probeFaces (" << probeFaces.size() << ")." << nl
            << "Number of patch probes and patch faces must match." << nl
            << exit(FatalError);
    }


    // Create and register the output heat-flux fields
    // Maybe number of probes is changed; thus, no check for the pointer
    qConvPtr_.reset
    (
        new scalarIOField
        (
            IOobject
            (
                "qConv",
                mesh().time().timeName(),
                mesh()
            ),
            szProbes_
        )
    );

    qRadPtr_.reset
    (
        new scalarIOField
        (
            IOobject
            (
                "qRad",
                mesh().time().timeName(),
                mesh()
            ),
            szProbes_
        )
    );


    // Read and store the common properties of the gauges
    Tgauge_ = dict.getScalar("Tgauge");
    if (Tgauge_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Gauge temperature cannot be negative." << nl
            << "Tgauge: " << Tgauge_
            << exit(FatalIOError);
    }

    a_ = dict.getOrDefault<scalar>("absorptivity", 1);
    e_ = dict.getOrDefault<scalar>("emissivity", 1);
    if (a_ < 0 || a_ > 1 || e_ < 0 || e_ > 1)
    {
        FatalIOErrorInFunction(dict)
            << "The range of absorptivity and emissivity can only be [0,1]."
            << "Absorptivity: " << a_ << nl
            << "Emissivity: " << e_
            << exit(FatalIOError);
    }


    // Read names of various fields; also early check if the fields exist
    Tname_ = dict.getOrDefault<word>("T", "T");
    const auto* Tptr = mesh().cfindObject<volScalarField>(Tname_);
    if (!Tptr)
    {
        FatalErrorInFunction
            << "Temperature field, T, is not available."
            << nl << "T: " << Tname_
            << exit(FatalError);
    }

    alphatName_ = dict.getOrDefault<word>("alphat", "alphat");
    const auto* alphatPtr = mesh().cfindObject<volScalarField>(alphatName_);
    if (!alphatPtr)
    {
        FatalErrorInFunction
            << "Turbulent thermal diffusivity field, alphat, is not available."
            << nl << "alphat: " << alphatName_
            << exit(FatalError);
    }


    // Compute the radiant-exitance term of the radiative heat flux
    // The term is calculated once per simulation unless Tgauge is changed
    Me_ = calcRadiantExitance();


    // Identify the faces that are actively sampled - assuming no duplication
    activeFaces_.resize_fill(probeFaces.size(), false);

    forAll(probeFaces, facei)
    {
        if (probeFaces[facei] != -1)
        {
            activeFaces_[facei] = true;
        }
    }


    // Identify the indices of the cells that contain the probed patch faces
    cells_ = identifyProbeCells();


    // Prepare the output files
    if (writeFile::canResetFile())
    {
        writeFile::resetFile(objName());
    }

    if (writeFile::canWriteHeader())
    {
        writeFileHeader(file());
        writtenHeader_ = true;
    }

    bool convective = dict.getOrDefault<bool>("convective", false);
    if (UPstream::master() && convective && !convectiveFilePtr_)
    {
        convectiveFilePtr_ = newFileAtStartTime("convective");
        writeFileHeader(convectiveFilePtr_());
    }

    bool radiative = dict.getOrDefault<bool>("radiative", false);
    if (UPstream::master() && radiative && !radiativeFilePtr_)
    {
        radiativeFilePtr_ = newFileAtStartTime("radiative");
        writeFileHeader(radiativeFilePtr_());
    }


    writeFields_ = dict.getOrDefault<bool>("writeFields", false);


    return true;
}


bool Foam::wallHeatFluxModels::gauge::execute()
{
    if (!(calcConvectiveHeatFlux() && calcRadiativeHeatFlux()))
    {
        return false;
    }
    return true;
}


bool Foam::wallHeatFluxModels::gauge::write()
{
    if (UPstream::master())
    {
        Ostream& os = file();

        os  << mesh().time().timeOutputValue();
        for (label pi = 0; pi < szProbes_; ++pi)
        {
            os  << ',' << (*qConvPtr_)[pi] + (*qRadPtr_)[pi];
        }
        os  << endl;
    }

    if (UPstream::master() && convectiveFilePtr_)
    {
        OFstream& os = convectiveFilePtr_.ref();

        os  << mesh().time().timeOutputValue();
        for (label pi = 0; pi < szProbes_; ++pi)
        {
            os  << ',' << (*qConvPtr_)[pi];
        }
        os  << endl;
    }

    if (UPstream::master() && radiativeFilePtr_)
    {
        OFstream& os = radiativeFilePtr_.ref();

        os  << mesh().time().timeOutputValue();
        for (label pi = 0; pi < szProbes_; ++pi)
        {
            os  << ',' << (*qRadPtr_)[pi];
        }
        os  << endl;
    }

    if (writeFields_)
    {
        if (qConvPtr_) qConvPtr_->write();
        if (qRadPtr_) qRadPtr_->write();
    }

    return true;
}


// ************************************************************************* //
