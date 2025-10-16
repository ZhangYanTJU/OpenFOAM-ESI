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

#include "FriedrichModel.H"
#include "cornerDetectionModel.H"
#include "processorFaPatch.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace filmSeparationModels
{
    defineTypeNameAndDebug(FriedrichModel, 0);
    addToRunTimeSelectionTable(filmSeparationModel, FriedrichModel, dictionary);


const Foam::Enum
<
    FriedrichModel::separationType
>
FriedrichModel::separationTypeNames
({
    { separationType::FULL, "full" },
    { separationType::PARTIAL , "partial" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

tmp<scalarField> FriedrichModel::Fratio() const
{
    const areaVectorField Up(film().Up());
    const areaVectorField& Uf = film().Uf();
    const areaScalarField& h = film().h();
    const areaScalarField& rho = film().rho();
    const areaScalarField& mu = film().mu();
    const areaScalarField& sigma = film().sigma();

    // Identify the faces where separation may occur
    const bitSet& separationFaces = cornerDetectorPtr_->getCornerFaces();

    // Calculate the corner angles corresponding to the separation faces
    const scalarList& separationAngles = cornerDetectorPtr_->getCornerAngles();


    // Initialize the force ratio
    auto tFratio = tmp<scalarField>::New(mesh().faces().size(), Zero);
    auto& Fratio = tFratio.ref();

    // Process internal faces
    forAll(separationFaces, i)
    {
        // Skip the routine if the face is not a candidate for separation
        if (!separationFaces[i]) continue;

        // Calculate the corner-angle trigonometric values
        const scalar sinAngle = std::sin(separationAngles[i]);
        const scalar cosAngle = std::cos(separationAngles[i]);

        // Reynolds number (FLW:Eq. 16)
        const scalar Re = h[i]*mag(Uf[i])*rho[i]/mu[i];

        // Weber number (FLW:Eq. 17)
        const vector Urel(Up[i] - Uf[i]);
        const scalar We = h[i]*rhop_*sqr(mag(Urel))/(2.0*sigma[i]);

        // Characteristic breakup length (FLW:Eq. 15)
        const scalar Lb =
            0.0388*Foam::sqrt(h[i])*Foam::pow(Re, 0.6)*Foam::pow(We, -0.5);

        // Force ratio - denominator (FLW:Eq. 20)
        const scalar den =
            sigma[i]*(sinAngle + 1.0) + rho[i]*magG_*h[i]*Lb*cosAngle;

        if (mag(den) > 0)
        {
            // Force ratio (FLW:Eq. 20)
            Fratio[i] = rho[i]*sqr(mag(Uf[i]))*h[i]*sinAngle/den;
        }
    }


    // Skip the rest of the routine if the simulation is a serial run
    if (!Pstream::parRun()) return tFratio;

    // Process processor faces
    const faBoundaryMesh& patches = mesh().boundary();

    for (const faPatch& fap : patches)
    {
        if (isA<processorFaPatch>(fap))
        {
            const label patchi = fap.index();
            const auto& edgeFaces = fap.edgeFaces();

            const auto& hp = h.boundaryField()[patchi];
            const auto& Ufp = Uf.boundaryField()[patchi];
            const auto& Upp = Up.boundaryField()[patchi];
            const auto& rhop = rho.boundaryField()[patchi];
            const auto& sigmap = sigma.boundaryField()[patchi];
            const auto& mup = mu.boundaryField()[patchi];

            forAll(hp, i)
            {
                // Skip the routine if the face is not a candidate for separation
                if (!separationFaces[i]) continue;

                const label faceO = edgeFaces[i];

                // Calculate the corner-angle trigonometric values
                const scalar sinAngle = std::sin(separationAngles[faceO]);
                const scalar cosAngle = std::cos(separationAngles[faceO]);

                // Reynolds number (FLW:Eq. 16)
                const scalar Re = hp[i]*mag(Ufp[i])*rhop[i]/mup[i];

                // Weber number (FLW:Eq. 17)
                const vector Urel(Upp[i] - Ufp[i]);
                const scalar We = hp[i]*rhop_*sqr(mag(Urel))/(2.0*sigmap[i]);

                // Characteristic breakup length (FLW:Eq. 15)
                const scalar Lb =
                    0.0388*Foam::sqrt(hp[i])
                   *Foam::pow(Re, 0.6)*Foam::pow(We, -0.5);

                // Force ratio - denominator (FLW:Eq. 20)
                const scalar den =
                    sigmap[i]*(sinAngle + 1.0)
                  + rhop[i]*magG_*hp[i]*Lb*cosAngle;

                if (mag(den) > 0)
                {
                    // Force ratio (FLW:Eq. 20)
                    Fratio[i] = rhop[i]*sqr(mag(Ufp[i]))*hp[i]*sinAngle/den;
                }
            }
        }
    }

    return tFratio;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

FriedrichModel::FriedrichModel
(
    const regionModels::areaSurfaceFilmModels::liquidFilmBase& film,
    const dictionary& dict
)
:
    filmSeparationModel(film, dict),
    separation_
    (
        separationTypeNames.getOrDefault
        (
            "separationType",
            dict,
            separationType::FULL
        )
    ),
    cornerDetectorPtr_(cornerDetectionModel::New(mesh(), film, dict)),
    rhop_(dict.getScalar("rhop")),
    magG_(mag(film.g().value())),
    C0_(dict.getOrDefault<scalar>("C0", 0.882)),
    C1_(dict.getOrDefault<scalar>("C1", -1.908)),
    C2_(dict.getOrDefault<scalar>("C2", 1.264))
{
    if (rhop_ < VSMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Primary-phase density, rhop: " << rhop_ << " must be non-zero."
            << abort(FatalIOError);
    }

    if (mag(C2_) < VSMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Empirical constant, C2 = " << C2_ << "cannot be zero."
            << abort(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

FriedrichModel::~FriedrichModel()
{}  // cornerDetectionModel was forward declared


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> FriedrichModel::separatedMassRatio() const
{
    cornerDetectorPtr_->detectCorners();

    tmp<scalarField> tFratio = Fratio();
    const auto& Fratio = tFratio.cref();

    // Initialize the mass ratio of film separation
    auto tseparated = tmp<scalarField>::New(mesh().faces().size(), Zero);
    auto& separated = tseparated.ref();


    switch (separation_)
    {
        case separationType::FULL:
        {
            forAll(Fratio, i)
            {
                if (Fratio[i] > 1)
                {
                    separated[i] = 1;
                }
            }
            break;
        }
        case separationType::PARTIAL:
        {
            forAll(Fratio, i)
            {
                if (Fratio[i] > 1)
                {
                    // (ZJD:Eq. 16)
                    separated[i] = C0_ + C1_*Foam::exp(-Fratio[i]/C2_);
                }
            }
            break;
        }
        default:
            break;  // This should not happen.
    }

    if (debug && mesh().time().writeTime())
    {
        {
            areaScalarField areaFratio
            (
                mesh().newIOobject("Fratio"),
                mesh(),
                dimensionedScalar(dimForce, Zero)
            );
            areaFratio.primitiveFieldRef() = Fratio;
            areaFratio.write();
        }

        {
            areaScalarField cornerAngles
            (
                mesh().newIOobject("cornerAngles"),
                mesh(),
                dimensionedScalar(dimless, Zero)
            );

            const bitSet& cornerFaces = cornerDetectorPtr_->getCornerFaces();
            const scalarList& angles = cornerDetectorPtr_->getCornerAngles();

            forAll(cornerFaces, i)
            {
                if (!cornerFaces[i]) continue;
                cornerAngles[i] = radToDeg(angles[i]);
            }
            cornerAngles.write();
        }
    }


    return tseparated;
}


/*
bool FriedrichModel::read(const dictionary& dict) const
{
    // Add the base-class reading later
    // Read the film separation model dictionary

    return true;
}
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace filmSeparationModels
} // End namespace Foam


// ************************************************************************* //

