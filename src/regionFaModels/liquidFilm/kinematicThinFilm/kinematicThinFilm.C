/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "kinematicThinFilm.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "gravityMeshObject.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kinematicThinFilm, 0);

addToRunTimeSelectionTable(liquidFilmBase, kinematicThinFilm, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool kinematicThinFilm::read(const dictionary& dict)
{
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kinematicThinFilm::kinematicThinFilm
(
    const word& modelType,
    const fvPatch& patch,
    const dictionary& dict
)
:
    liquidFilmModel(modelType, patch, dict)
{
    init();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

kinematicThinFilm::~kinematicThinFilm()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void kinematicThinFilm::init()
{
}


void kinematicThinFilm::preEvolveRegion()
{
    // Update reads
    liquidFilmModel::read(coeffs());

    // Update mass exchange sources
    liquidFilmModel::preEvolveRegion();

    // gas pressure map from primary region
    ppf_ = pg();
}


void kinematicThinFilm::evolveRegion()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    const uniformDimensionedVectorField& g=meshObjects::gravity::New(time());

    const areaVectorField gs(g - gn_*regionMesh().faceAreaNormals());

    edgeScalarField phi2s("phi2s", fac::interpolate(h_)*phif_);

    // Solve continuity for h
    for (int oCorr=1; oCorr<=nOuterCorr_; oCorr++)
    {
        faVectorMatrix UsEqn
        (
            fam::ddt(h_, Uf_)
          + fam::div(phi2s, Uf_)
          ==
            gs*h_
          + turbulence_->Su(Uf_)
          + faOptions()(h_, Uf_, sqr(dimVelocity))
          + rhoSp_*Uf_
          + USp_
        );

        faOptions().constrain(UsEqn);

        UsEqn.relax();

        if (momentumPredictor_)
        {
            solve(UsEqn == -fac::grad(pf_*h_)/rho_ + pf_*fac::grad(h_)/rho_);
        }

        for (int corr=1; corr<=nCorr_; corr++)
        {
            areaScalarField UsA(UsEqn.A());

            Uf_ = UsEqn.H()/UsA;
            Uf_.correctBoundaryConditions();

            phif_ =
                (fac::interpolate(Uf_) & regionMesh().Le())
                - fac::interpolate(1.0/(rho_*UsA))
                * fac::lnGrad(pf_*h_)*regionMesh().magLe()
                + fac::interpolate(pf_/(rho_*UsA))
                * fac::lnGrad(h_)*regionMesh().magLe();

            for (int nonOrth=0; nonOrth<=nNonOrthCorr_; nonOrth++)
            {
                faScalarMatrix hEqn
                (
                    fam::ddt(h_)
                  + fam::div(phif_, h_)
                 ==
                    faOptions()(rho_, h_, dimVelocity)
                  + rhoSp_
                );

                faOptions().constrain(hEqn);

                hEqn.relax();
                hEqn.solve();

                if (nonOrth == nNonOrthCorr_)
                {
                    phi2s = hEqn.flux();
                }
            }

            // Bound h_
            h_.primitiveFieldRef() = max
            (
                max
                (
                    h_.primitiveField(),
                    fac::average(max(h_, h0_))().primitiveField()
                    *pos(h0_.value() - h_.primitiveField())
                ),
                h0_.value()
            );

            pf_ = rho_*gn_*h_ - sigma_*fac::laplacian(h_) + ppf_ + pnSp_;
            pf_.correctBoundaryConditions();

            Uf_ -= (1.0/(rho_*UsA))*fac::grad(pf_*h_)
                 - (pf_/(rho_*UsA))*fac::grad(h_);
            Uf_.correctBoundaryConditions();
        }
    }

    // Update deltaRho_ with new delta_
    //hRho_ == h_*rho_;
}

void kinematicThinFilm::postEvolveRegion()
{
    // Reset sources
    liquidFilmModel::postEvolveRegion();

    // Correct thermo
    correctThermoFields();

    // Correct turbulence
    turbulence_->correct();

}


void kinematicThinFilm::info()
{
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
