/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024-2025 OpenCFD Ltd.
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

#include "bladeForces.H"
#include "fvcGrad.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "cartesianCS.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"
#include "polySurfaceFields.H"
// #include "foamVtkSurfaceWriter.H"
// #include "surfaceWriter.H"
#include "vtkSurfaceWriter.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(bladeForces, 0);
    addToRunTimeSelectionTable(functionObject, bladeForces, dictionary);
}
}


// Somewhat experimental...
#undef  Foam_bladeForces_CalculateChord

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

#ifdef Foam_bladeForces_CalculateChord

// Reduce min/max values - MUST have identical length on all ranks!
// Uses en-mass reduction
void reduceMinMax
(
    List<scalarMinMax>& limits,
    List<scalar>& work,  // work array
    const int tag = UPstream::msgType(),
    const label comm = UPstream::worldComm
)
{
    const label len = limits.size();

    work.resize_nocopy(len);

    for (label i = 0; i < len; ++i)
    {
        work[i] = limits[i].min();
    }

    Foam::reduce(work.data(), work.size(), minOp<scalar>(), tag, comm);

    for (label i = 0; i < len; ++i)
    {
        limits[i].min() = work[i];
        work[i] = limits[i].max();
    }

    Foam::reduce(work.data(), work.size(), maxOp<scalar>(), tag, comm);

    for (label i = 0; i < len; ++i)
    {
        limits[i].max() = work[i];
    }
}

#endif  /* Foam_bladeForces_CalculateChord */

} // End namespace Foam


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::bladeForces::initialise()
{
    if (initialised_)
    {
        return;
    }

    const bool isCompressible = (rhoName_ != "rhoInf");

    if
    (
        !foundObject<volVectorField>(UName_)
     || !foundObject<volScalarField>(pName_)
     || (isCompressible && !foundObject<volScalarField>(rhoName_))
    )
    {
        FatalErrorInFunction
            << "Could not find velocity '" << UName_
            << "' or pressure '" << pName_ << '\'';

        if (isCompressible)
        {
            FatalError
                << " or density '" << rhoName_ << '\'';
        }

        FatalError
            << " in database" << endl
            << exit(FatalError);
    }

    initialised_ = true;
}


Foam::tmp<Foam::symmTensorField>
Foam::functionObjects::bladeForces::devRhoReff
(
    const tensorField& gradUp,
    const label patchi
) const
{
    if
    (
        typedef incompressible::turbulenceModel icoTurbModel;
        const auto* turbp
      = cfindObject<icoTurbModel>(icoTurbModel::propertiesName)
    )
    {
        const auto& turb = *turbp;

        return -rho(patchi)*turb.nuEff(patchi)*devTwoSymm(gradUp);
    }

    if
    (
        typedef compressible::turbulenceModel cmpTurbModel;
        const auto* turbp
      = cfindObject<cmpTurbModel>(cmpTurbModel::propertiesName)
    )
    {
        const auto& turb = *turbp;

        return -turb.muEff(patchi)*devTwoSymm(gradUp);
    }

    if
    (
        const auto* thermop
      = cfindObject<fluidThermo>(fluidThermo::dictName)
    )
    {
        const auto& thermo = *thermop;

        return -thermo.mu(patchi)*devTwoSymm(gradUp);
    }

    if
    (
        const auto* props
      = cfindObject<transportModel>("transportProperties")
    )
    {
        const auto& laminarT = *props;

        return -rho(patchi)*laminarT.nu(patchi)*devTwoSymm(gradUp);
    }

    if
    (
        const auto* props
      = cfindObject<dictionary>("transportProperties")
    )
    {
        const dimensionedScalar nu("nu", dimViscosity, *props);

        return -rho(patchi)*nu.value()*devTwoSymm(gradUp);
    }
    else
    {
        // Failed to resolve any model
        FatalErrorInFunction
            << "No valid model for viscous stress calculation"
            << exit(FatalError);

        return nullptr;
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::bladeForces::mu() const
{
    if
    (
        const auto* thermop
      = cfindObject<fluidThermo>(basicThermo::dictName)
    )
    {
        const auto& thermo = *thermop;

        return thermo.mu();
    }

    if
    (
        const auto* props
      = cfindObject<transportModel>("transportProperties")
    )
    {
        const auto& laminarT = *props;

        return rho()*laminarT.nu();
    }

    if
    (
        const auto* props
      = cfindObject<dictionary>("transportProperties")
    )
    {
        const dimensionedScalar nu("nu", dimViscosity, *props);

        return rho()*nu;
    }
    else
    {
        // Failed to resolve any model
        FatalErrorInFunction
            << "No valid model for dynamic viscosity calculation"
            << exit(FatalError);

        return nullptr;
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::bladeForces::rho() const
{
    if (rhoName_ == "rhoInf")
    {
        return volScalarField::New
        (
            "rho",
            mesh_,
            dimensionedScalar(dimDensity, rhoRef_)
        );
    }

    return lookupObject<volScalarField>(rhoName_);
}


Foam::tmp<Foam::scalarField>
Foam::functionObjects::bladeForces::rho(const label patchi) const
{
    if (rhoName_ == "rhoInf")
    {
        return tmp<scalarField>::New
        (
            mesh_.boundary()[patchi].size(),
            rhoRef_
        );
    }

    const auto& rho = lookupObject<volScalarField>(rhoName_);
    return rho.boundaryField()[patchi];
}


Foam::scalar Foam::functionObjects::bladeForces::rho
(
    const volScalarField& p
) const
{
    if (p.dimensions() == dimPressure)
    {
        return 1;
    }

    if (rhoName_ != "rhoInf")
    {
        FatalErrorInFunction
            << "Dynamic pressure is expected but kinematic is provided."
            << exit(FatalError);
    }

    return rhoRef_;
}


void Foam::functionObjects::bladeForces::createIntegratedDataFiles()
{
    if (!forceFilePtr_)
    {
        forceFilePtr_ = newFileAtStartTime("force");
        auto& os = forceFilePtr_();

        writeHeader(os, "Force");
        writeHeader(os, "");
        writeCommented(os, "Time");
        writeTabbed(os, "thrust");
        writeTabbed(os, "drag");
        writeTabbed(os, "torque");
        writeTabbed(os, "Cd");
        writeTabbed(os, "Cl");
        writeTabbed(os, "Cp");
        os  << endl;
    }

    const auto writeCoeffsHeader =
        [&](auto& os, const auto& longName, const auto& shortName)
        {
            writeHeader(os, longName);
            writeHeader(os, "");
            writeCommented(os, "number of radial bins:");
            os << ' ' << nRadialDiv_ << nl;

            // General formatting (not fixed/scientific) for the bin limits:
            const auto oldFlags = os.flags();
            os.unsetf(std::ios_base::floatfield);

            writeCommented(os, "bin lower limits:");
            for (label bin = 0; bin < nRadialDiv_; ++bin)
            {
                os << ' ' << float(bin)/nRadialDiv_;
            }
            os << nl;
            writeCommented(os, "bin upper limits:");
            for (label bin = 1; bin <= nRadialDiv_; ++bin)
            {
                os << ' ' << float(bin)/nRadialDiv_;
            }
            os  << nl;
            writeHeader(os, "");
            writeCommented(os, "Time");
            writeTabbed(os, shortName);
            os  << endl;

            os.setf(oldFlags);
        };


    if (!CdFilePtr_)
    {
        CdFilePtr_ = newFileAtStartTime("C_d");
        auto& os = CdFilePtr_();

        writeCoeffsHeader(os, "Drag coefficient", "Cd(per bin)");
    }

    if (!ClFilePtr_)
    {
        ClFilePtr_ = newFileAtStartTime("C_l");
        auto& os = ClFilePtr_();

        writeCoeffsHeader(os, "Lift coefficient", "Cl(per bin)");
    }

    if (!CpFilePtr_)
    {
        CpFilePtr_ = newFileAtStartTime("C_p");
        auto& os = CpFilePtr_();

        writeCoeffsHeader(os, "Pressure coefficient", "Cp(per bin)");
    }
}


void Foam::functionObjects::bladeForces::writeIntegratedDataFiles()
{
    {
        auto& os = forceFilePtr_();

        writeCurrentTime(os);

        writeValue(os, (useMagThrust_ ? Foam::mag(sumThrust_) : sumThrust_));
        writeValue(os, (useMagDrag_ ? Foam::mag(sumDrag_) : sumDrag_));
        writeValue(os, sumTorque_);
        writeValue(os, totalCd_);
        writeValue(os, totalCl_);
        writeValue(os, totalCp_);
        os  << endl;
    }

    {
        auto& os = CdFilePtr_();

        writeCurrentTime(os);

        for (const auto& val : bandCd_)
        {
            writeValue(os, val);
        }
        os  << endl;
    }

    {
        auto& os = ClFilePtr_();

        writeCurrentTime(os);

        for (const auto& val : bandCl_)
        {
            writeValue(os, val);
        }
        os  << endl;
    }

    {
        auto& os = CpFilePtr_();

        writeCurrentTime(os);

        for (const auto& val : bandCp_)
        {
            writeValue(os, val);
        }
        os  << endl;
    }
}


void Foam::functionObjects::bladeForces::logIntegratedData
(
    const string& descriptor,
    const vector& pres,
    const vector& vis,
    const vector& internal
) const
{
//    if (!log)
//    {
//        return;
//    }
//
//    Log << "    Sum of " << descriptor.c_str() << nl
//        << "        Total    : " << (pres + vis + internal) << nl
//        << "        Pressure : " << pres << nl
//        << "        Viscous  : " << vis << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::bladeForces::bladeForces
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    bool readFields
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name),

    sumThrust_(0),
    sumDrag_(0),
    sumTorque_(0),

    totalCd_(0),
    totalCl_(0),
    totalCp_(0),

    refRadius_(1),
    nRadialDiv_(10),
    revPerSec_(1),

    Uref_(1),
    pRef_(0),
    rhoRef_(1),
    pName_("p"),
    UName_("U"),
    rhoName_("rho"),
    outputName_("blade"),
    writeCounter_(0),
    fieldsInterval_(0),
    initialised_(false),
    writeFields_(false),
    useGeometricVelocity_(false),
    useNearCellValue_(false),
    lefthand_(false),
    useMagThrust_(false),
    useMagDrag_(false)
{
    if (readFields)
    {
        read(dict);
        Log << endl;
    }
}


Foam::functionObjects::bladeForces::bladeForces
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    bool readFields
)
:
    fvMeshFunctionObject(name, obr, dict),
    writeFile(mesh_, name),

    sumThrust_(0),
    sumDrag_(0),
    sumTorque_(0),

    totalCd_(0),
    totalCl_(0),
    totalCp_(0),

    refRadius_(1),
    nRadialDiv_(10),
    revPerSec_(1),

    Uref_(1),
    pRef_(0),
    rhoRef_(1),
    pName_("p"),
    UName_("U"),
    rhoName_("rho"),
    outputName_("blade"),
    writeCounter_(0),
    fieldsInterval_(0),
    initialised_(false),
    writeFields_(false),
    useGeometricVelocity_(false),
    useNearCellValue_(false),
    lefthand_(false),
    useMagThrust_(false),
    useMagDrag_(false)
{
    if (readFields)
    {
        read(dict);
        Log << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::bladeForces::read(const dictionary& dict)
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    if (!fvMeshFunctionObject::read(dict) || !writeFile::read(dict))
    {
        return false;
    }

    initialised_ = false;

    Info<< type() << ' ' << name() << ':' << nl;

    // Can also use pbm.indices(), but no warnings...
    patchIDs_ = pbm.patchSet(dict.get<wordRes>("patches")).sortedToc();

    refRadius_ = dict.getOrDefault<scalar>("radius", 1);
    nRadialDiv_ = dict.getOrDefault<label>("nRadial", 10);
    revPerSec_ = dict.get<scalar>("n");

    cylCoord_ = coordSystem::cylindrical
    (
        dict.get<point>("origin"),
        dict.get<vector>("axis")
    );

    {
        // Optional field name entries
        if (dict.readIfPresent("p", pName_))
        {
            Info<< "    p: " << pName_ << nl;
        }
        if (dict.readIfPresent("U", UName_))
        {
            Info<< "    U: " << UName_ << nl;
        }
        if (dict.readIfPresent("rho", rhoName_))
        {
            Info<< "    rho: " << rhoName_ << nl;
        }

        // Reference axial velocity, 1 by default
        Uref_ = dict.get<scalar>("Uref");
        {
            Info<< "    Reference axial velocity (Uref) set to "
                << Uref_ << nl;
        }

        // Reference density (1 by default) for incompressible calculations
        rhoRef_ = 1;
        if
        (
            rhoName_ == "rhoInf"
         && dict.readIfPresent("rhoInf", rhoRef_)
        )
        {
            Info<< "    Freestream density (rhoInf) set to "
                << rhoRef_ << nl;
        }

        // Reference pressure, 0 by default
        if (dict.readIfPresent("pRef", pRef_))
        {
            Info<< "    Reference pressure (pRef) set to " << pRef_ << nl;
        }
    }

    // Name for VTP output and local surfaces
    dict.readIfPresent("outputName", outputName_);

    useNearCellValue_ = false;
    useGeometricVelocity_ = false;

    if
    (
        dict.readIfPresent("nearCellValue", useNearCellValue_)
     && useNearCellValue_
    )
    {
        Info<< "    Blade wall velocity from near cell values" << nl;
    }
    else if
    (
        dict.readIfPresent("geometricVelocity", useGeometricVelocity_)
     && useGeometricVelocity_
    )
    {
        Info<< "    Blade wall velocity from geometry" << nl;
    }


    writeCounter_ = 0;
    fieldsInterval_ = 0;
    dict.readIfPresent("fieldsInterval", fieldsInterval_);

    if
    (
        writeFields_ = false;
        dict.readIfPresent("writeFields", writeFields_) && writeFields_
    )
    {
        Info<< "    Fields will be written ";

        if (fieldsInterval_ <= 1)
        {
            Info<< "each write interval" << nl;
        }
        else
        {
            Info<< "every " << fieldsInterval_ << " write intervals" << nl;
        }
    }

    if
    (
        lefthand_ = false;
        dict.readIfPresent("lefthand", lefthand_) && lefthand_
    )
    {
        Info<< "    Left-hand blade = true" << nl;
    }

    if
    (
        useMagThrust_ = false;
        dict.readIfPresent("mag.thrust", useMagThrust_) && useMagThrust_
    )
    {
        Info<< "    Ignore sign for thrust" << nl;
    }

    if
    (
        useMagDrag_ = false;
        dict.readIfPresent("mag.drag", useMagDrag_) && useMagDrag_
    )
    {
        Info<< "    Ignore sign for drag" << nl;
    }

    Info<< "    Output surface name: " << outputName_ << endl;

    return true;
}


// Internal fields on the "blade" (outputName) surface:
//
// - "U", "p" :
// - "force" :
// - "Cd" :
// - "Cl" :
// - "Cp" :
// - "radius" : The dimensionless radius
// - "radial-bins" : The radial bin ID for each selected patch face
// - "patch-ids" : The patch IDs. Can use to apply a threshold filter

void Foam::functionObjects::bladeForces::calculate()
{
    initialise();

    sumThrust_ = 0;
    sumDrag_ = 0;
    sumTorque_ = 0;

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    // Define (or redefine) the blade surface, saved to storedObjects()
    // to allow further post-processing (eg, with surfaceFieldValue)

    auto& surf = polySurface::New(scopedName(outputName_), storedObjects());

    // Pout<< "patchIDs : " << patchIDs_ << nl;

    // Offsets within compact addressing
    labelList compactOffsets(patchIDs_.size()+1);

    {
        // Offsets for number of faces
        compactOffsets[0] = 0;

        forAll(patchIDs_, i)
        {
            const polyPatch& pp = pbm[patchIDs_[i]];

            compactOffsets[i+1] = compactOffsets[i] + pp.size();
        }

        labelList meshFaceIds(compactOffsets.back());

        forAll(patchIDs_, i)
        {
            const polyPatch& pp = pbm[patchIDs_[i]];

            auto slice = meshFaceIds.slice
            (
                compactOffsets[i],
                (compactOffsets[i+1] - compactOffsets[i])
            );

            Foam::identity(slice, pp.start());
        }

        // Reset geometry (and fields)

        indirectPrimitivePatch allPatches
        (
            IndirectList<face>(mesh_.faces(), std::move(meshFaceIds)),
            mesh_.points()
        );

        surf.clearFields();
        surf.copySurface(allPatches.localPoints(), allPatches.localFaces());
    }

    // Pout<< "surf : " << surf.nFaces() << nl;

    // Assign the correct pressure dimensions later on
    auto& force = surf.newField<vector>("force", dimForce);
    auto& p_surf = surf.newField<scalar>("p", dimless);
    auto& U_surf = surf.newField<vector>("U", dimVelocity);
    auto& Cd = surf.newField<scalar>("Cd", dimless);
    auto& Cl = surf.newField<scalar>("Cl", dimless);
    auto& Cp = surf.newField<scalar>("Cp", dimless);

    // The dimensionless radius
    auto& radius = surf.newField<scalar>("radius", dimless);

    // The radial bin for each selected patch face
    auto& radialBins = surf.newField<label>("radial-bins", dimless);

    // The patch id for selected patch face
    auto& patchBins = surf.newField<label>("patch-ids", dimless);


    scalarField faceAreas(surf.nFaces(), Foam::zero{});

    // Or local only
    /// The radial bin for each selected patch face
    /// labelList radialBins(surf.nFaces(), Foam::zero{});

    // scalarList thrust(nRadialDiv_, Foam::zero{});
    // scalarList torque(nRadialDiv_, Foam::zero{});

    // Coarse estimation of the chord length
    // - find min/max angle for each radial bin,
    // - apply with the centre-line radius
    //
    // Warnings:
    // * this length is effectively the projected one (without depth).
    // * angle determination may be overestimated if the faces on the upper
    //   radius are larger than a simple sweep.
    // * angle will be underestimated by 1/2 face size if there are no
    //   leading/trailing surfaces included
    // * does not work really at all for multiple blades
    // * needs special (ad hoc) handling of -pi/+pi quadrant crossing

    scalarList bandArea(nRadialDiv_, Foam::zero{});

    #ifdef Foam_bladeForces_CalculateChord
    scalarList chordLength(nRadialDiv_, Foam::zero{});
    #endif

    bandCd_.resize_nocopy(nRadialDiv_);
    bandCd_ = Foam::zero{};

    bandCl_.resize_nocopy(nRadialDiv_);
    bandCl_ = Foam::zero{};

    bandCp_.resize_nocopy(nRadialDiv_);
    bandCp_ = Foam::zero{};


    // The tangential drag direction
    const int dragDirection = (lefthand_ ? -1 : +1);

    // From rev/sec -> 1/s
    const scalar omega = (constant::mathematical::twoPi * revPerSec_);

    // Step 0:
    // - extract geometric information
    // - populate geometric velocity field (optional)

    {
        const auto& Cb = mesh_.C().boundaryField();
        const auto& Sfb = mesh_.Sf().boundaryField();

        // The min/max limits for calculating chord length
        #ifdef Foam_bladeForces_CalculateChord
        List<scalarMinMax> chordLimitsAngle(nRadialDiv_);
        List<scalarMinMax> chordLimitsDepth(nRadialDiv_);
        #endif

        // Classify relative radius to radial bin number
        const auto whichRadialBin = [=](scalar relRadius) -> label
        {
            const label bin = std::floor(nRadialDiv_ * relRadius);

            return
            (
                (bin <= 0)
              ? 0
              : (bin >= nRadialDiv_)
              ? (nRadialDiv_ - 1)
              : bin
            );
        };

        forAll(patchIDs_, zonei)
        {
            const label patchi = patchIDs_[zonei];

            // Sub-slice within sampled surface fields
            const labelRange range
            (
                compactOffsets[zonei],
                compactOffsets[zonei+1] - compactOffsets[zonei]
            );

            const label nPatchFaces = range.size();
            // Same as  nPatchFaces = pbm[patchi].size();

            if (!nPatchFaces)
            {
                continue;
            }

            // The patch id
            patchBins.slice(range) = patchi;


            // Local cylindrical (r-theta-z) positions
            const pointField localPosition
            (
                cylCoord_.localPosition(Cb[patchi])
            );

            auto binsb = radialBins.slice(range);
            auto radiusb = radius.slice(range);
            auto areab = faceAreas.slice(range);

            // Optionally calculate blade speed from geometry
            if (useGeometricVelocity_)
            {
                auto Ub = U_surf.field().slice(range);
                const auto& patchCentres = Cb[patchi];
                const auto& patchSf = Sfb[patchi];

                forAll(localPosition, facei)
                {
                    Ub[facei] =
                    (
                        cylCoord_.transform
                        (
                            patchCentres[facei],
                            vector(0, localPosition[facei].x()*omega, 0)
                        )
                    );

                    // Eliminate normal components
                    Ub[facei].removeCollinear(normalised(patchSf[facei]));
                }
            }

            forAll(localPosition, facei)
            {
                const auto& cylPoint = localPosition[facei];

                radiusb[facei] = (cylPoint.x()/refRadius_);

                const label binId = whichRadialBin(radiusb[facei]);

                binsb[facei] = binId;

                areab[facei] = Foam::mag(Sfb[patchi][facei]);
                bandArea[binId] += areab[facei];

                #ifdef Foam_bladeForces_CalculateChord
                chordLimitsAngle[binId].add(cylPoint.y());
                chordLimitsDepth[binId].add(cylPoint.z());
                #endif
            }
        }

        // Reductions - perform en-mass (identical length on all procs)

        Foam::reduce
        (
            bandArea.data(),
            bandArea.size(),
            sumOp<scalar>(),
            UPstream::msgType(),
            UPstream::worldComm
        );

        #ifdef Foam_bladeForces_CalculateChord
        {
            scalarList values;  // work array
            reduceMinMax(chordLimitsAngle, values);
            reduceMinMax(chordLimitsDepth, values);

            // Have -pi/+pi range. Crossing between the (-x/-y) and (-x/+y)
            // quadrants will always (-pi/+pi) limits and thus does not
            // correspond to the blade edges - ie, useless for calculating
            // chord length.
            //
            // Ad hoc fix is to detect if such a crossing occurs
            // and rotate the quadrants accordingly.

            bool crossesQuadrant(false);

            for (const auto& angles : chordLimitsAngle)
            {
                if (angles.min() < -3.14 && angles.max() > 3.14)
                {
                    crossesQuadrant = true;
                    break;
                }
            }

            if (crossesQuadrant)
            {
                // Blade crosses between (-x/-y) and (-x/+y) quadrants

                chordLimitsAngle = scalarMinMax();

                for (const label patchi : patchIDs_)
                {
                    // Local cylindrical (r-theta-z) positions
                    const pointField localPosition
                    (
                        cylCoord_.localPosition(Cb[patchi])
                    );

                    for (const auto& cylPoint : localPosition)
                    {
                        const label binId =
                            whichRadialBin(cylPoint.x()/refRadius_);

                        scalar theta = cylPoint.y();

                        if (theta < 0)
                        {
                            // Rotate (-x/-y) quadrant to (+x/+y)
                            theta += constant::mathematical::pi;
                        }
                        else
                        {
                            // Rotate (-x/+y) quadrant to (+x/-y)
                            theta -= constant::mathematical::pi;
                        }

                        chordLimitsAngle[binId].add(theta);
                    }
                }

                reduceMinMax(chordLimitsAngle, values);
            }

            forAll(chordLength, i)
            {
                if (chordLimitsAngle[i].empty() || chordLimitsDepth[i].empty())
                {
                    chordLength[i] = 0;
                    continue;
                }

                const scalar angle = chordLimitsAngle[i].span();
                const scalar depth = chordLimitsDepth[i].span();

                chordLength[i] = Foam::hypot
                (
                    // Use the mid-radius
                    (((0.5 + i)*refRadius_)/nRadialDiv_)*angle,
                    depth
                );
            }
        }
        #endif
    }

    // Now have per-band:
    // - bandArea (reduced)
    //
    // Now have per-face:
    // - bin, radius, area

    {
        // Pressure (volume, boundary)
        const auto& p = lookupObject<volScalarField>(pName_);
        const auto& pb = p.boundaryField();

        // Assign the correct pressure dimensions
        p_surf.dimensions().reset(p.dimensions());

        // Face centres, face area normals (boundary)
        const auto& Cb = mesh_.C().boundaryField();
        const auto& Sfb = mesh_.Sf().boundaryField();

        // Velocity (volume, boundary)
        const auto& U = lookupObject<volVectorField>(UName_);
        const auto& Up = U.boundaryField();

        const volTensorField gradU(fvc::grad(U));
        const auto& gradUb = gradU.boundaryField();

        // Scale pRef by density for incompressible simulations
        const scalar rhoRef = rho(p);
        const scalar pRef = pRef_/rhoRef;

        forAll(patchIDs_, zonei)
        {
            const label patchi = patchIDs_[zonei];

            // Sub-slice within sampled surface fields
            const labelRange range
            (
                compactOffsets[zonei],
                compactOffsets[zonei+1] - compactOffsets[zonei]
            );

            const label nPatchFaces = range.size();
            // Same as  nPatchFaces = pbm[patchi].size();

            if (!nPatchFaces)
            {
                continue;
            }

            const auto binsb = radialBins.slice(range);
            const auto radiusb = radius.slice(range);
            const auto areab = faceAreas.slice(range);

            // Density (on patch)
            const tmp<scalarField> trhop = rho(patchi);
            const auto& rhop = trhop();

            // Pressure (on patch)
            const auto& pressp = pb[patchi];

            // Add pressure to the surface output
            p_surf.field().slice(range) = pressp;

            // Add velocity (as near cell velocity) to the surface output
            if (false && (useNearCellValue_ && !useGeometricVelocity_))
            {
                auto Ub = U_surf.field().slice(range);

                // From geometry ...
                // if (prescribedVelocity_)
                //
                // calculate from geometry...

                if (useNearCellValue_)
                {
                    // Near cell velocity
                    const auto tUcellb = Up[patchi].patchInternalField();

                    if (Ub.size() == tUcellb().size())
                    {
                        Ub = tUcellb();
                    }
                }
                else if (!useGeometricVelocity_)
                {
                    // Patch values
                    if (Ub.size() == Up[patchi].size())
                    {
                        Ub = Up[patchi];
                    }
                }
            }

            auto Cdb = Cd.field().slice(range);
            auto Clb = Cl.field().slice(range);
            auto Cpb = Cp.field().slice(range);


            // Forces : global reference frame
            const vectorField globalForce
            (
                // Pressure forces
                (rhoRef*Sfb[patchi]*(pb[patchi] - pRef))

                // Viscous forces
              + (Sfb[patchi] & devRhoReff(gradUb[patchi], patchi))
            );

            // Forces : locally oriented (r-theta-z local vectors)
            const vectorField localForce
            (
                cylCoord_.invTransform(Cb[patchi], globalForce)
            );


            // Globally oriented
            force.field().slice(range) = globalForce;


            for (label facei = 0; facei < nPatchFaces; ++facei)
            {
                const scalar r = (radiusb[facei] * refRadius_);
                const label binId = binsb[facei];

                // thrust = force(axial)
                const scalar thrustOnFace = localForce[facei].z();

                // drag = force(tangential)
                const scalar dragOnFace = (dragDirection*localForce[facei].y());

                // torque = radius * force(tangential)
                const scalar torqueOnFace = (r*dragOnFace);

                sumDrag_   += dragOnFace;
                sumThrust_ += thrustOnFace;
                sumTorque_ += torqueOnFace;

                // Pressure coefficient:
                //     Cp = (p - pInf) / (1/2 * rho * Vr^2)
                //     where Vr^2 = Vx^2 + (2 pi*n * r)^2

                const scalar Urel2 = (sqr(Uref_) + sqr(omega*r));

                Cpb[facei] =
                (
                    (pressp[facei] - pRef) / (0.5*rhop[facei]*Urel2)
                );

                // Drag coefficient:
                //     Cd = (Drag / (1/2 * rho * V^2 * A))
                Cdb[facei] =
                (
                    dragOnFace / (0.5*rhop[facei]*Urel2*bandArea[binId])
                );

                // Lift coefficient:
                //     Cl = (Lift / (1/2 * rho * V^2 * A))
                Clb[facei] =
                (
                    thrustOnFace / (0.5*rhop[facei]*Urel2*bandArea[binId])
                );

                bandCd_[binId] += Cdb[facei]*areab[facei];
                bandCl_[binId] += Clb[facei]*areab[facei];
                bandCp_[binId] += Cpb[facei]*areab[facei];
            }
        }
    }

    // Global reduction followed by sum-area weighting
    // return the total area-weighted value.
    // Note: the bandArea has already been reduced above.
    const auto reduce_normalize = [&](UList<scalar>& bandField) -> scalar
    {
        Foam::reduce
        (
            bandField.data(),
            bandField.size(),
            sumOp<scalar>(),
            UPstream::msgType(),
            UPstream::worldComm
        );

        scalar totalWeighted(0);
        scalar totalArea(0);

        forAll(bandField, binId)
        {
            totalWeighted += bandField[binId];
            totalArea += bandArea[binId];

            bandField[binId] /= (bandArea[binId] + VSMALL);
        }

        return (totalWeighted / (totalArea + VSMALL));
    };

    totalCd_ = reduce_normalize(bandCd_);
    totalCl_ = reduce_normalize(bandCl_);
    totalCp_ = reduce_normalize(bandCp_);

    Foam::reduce(sumThrust_, sumOp<scalar>());
    Foam::reduce(sumDrag_, sumOp<scalar>());
    Foam::reduce(sumTorque_, sumOp<scalar>());
}


void Foam::functionObjects::bladeForces::writeSurface() const
{
    const auto* surfptr =
        storedObjects().cfindObject<polySurface>(scopedName(outputName_));

    if (!surfptr)
    {
        Info<< "missing surface " << scopedName(outputName_)
            << " in " << storedObjects().sortedNames() << nl;
        return;
    }
    const auto& surf = *surfptr;

    // Output may need more work...
    {
        surfaceWriters::vtkWriter writer;

        // Use outputDir/TIME/surface-name
        writer.useTimeDir(true);

        // Time-aware, with time spliced into the output path
        writer.beginTime(time_);

        writer.open
        (
            surf.points(),
            surf.faces(),
            (baseFileDir() / name() / "surfaces" / outputName_)
            // parallel = true
        );

        writer.nFields(surf.nFaceData());  // Legacy VTK

        #undef  doLocalCode
        #define doLocalCode(Type)                                             \
        {                                                                     \
            for (const auto& fld : surf.csorted<PolyFaceField<Type>>())       \
            {                                                                 \
                writer.write(fld.name(), fld);                                \
            }                                                                 \
        }

        doLocalCode(label);
        doLocalCode(scalar);
        doLocalCode(vector);
        #undef doLocalCode

        writer.endTime();
    }
}


bool Foam::functionObjects::bladeForces::execute()
{
    calculate();

    Log << type() << ' ' << name() << " :" << nl;

    Log << "    Thrust : "
        << (useMagThrust_ ? Foam::mag(sumThrust_) : sumThrust_) << nl
        << "    Drag   : "
        << (useMagDrag_ ? Foam::mag(sumDrag_) : sumDrag_) << nl
        << "    Torque : " << sumTorque_ << nl
        << "    Cd     : " << totalCd_ << nl
        << "    Cl     : " << totalCl_ << nl
        << "    Cp     : " << totalCp_ << nl;

    setResult("thrust", sumThrust_);
    setResult("drag", sumDrag_);
    setResult("torque", sumTorque_);

    return true;
}


bool Foam::functionObjects::bladeForces::write()
{
    bool hasOutput(false);

    ++writeCounter_;

    if (writeToFile())
    {
        hasOutput = true;
        Log << type() << ' ' << name()
            << " : writing force and coefficient files" << nl;

        createIntegratedDataFiles();
        writeIntegratedDataFiles();
    }

    if
    (
        writeFields_ &&
        ((fieldsInterval_ <= 1) || !(writeCounter_ % fieldsInterval_))
    )
    {
        hasOutput = true;
        Log << type() << ' ' << name()
            << " : writing surface and fields -> "
            << scopedName(outputName_) << nl;

        writeSurface();
    }

    if (hasOutput)
    {
        Log << endl;
    }

    return true;
}


// ************************************************************************* //
