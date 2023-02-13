/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "chargeCloud.H"

#include "makeReactingParcelCloudFunctionObjects.H"

// Kinematic
#include "makeThermoParcelForces.H" // thermo variant
#include "makeParcelDispersionModels.H"
#include "makeChargeParcelInjectionModels.H"
#include "makeParcelPatchInteractionModels.H"
#include "makeReactingMultiphaseParcelStochasticCollisionModels.H" // MP variant
#include "makeReactingParcelSurfaceFilmModels.H" // Reacting variant

// Thermodynamic
#include "makeParcelHeatTransferModels.H"

// Reacting
#include "makeReactingMultiphaseParcelCompositionModels.H" // MP Variant
#include "makeReactingParcelPhaseChangeModels.H"

// Reacting multiphase
#include "makeReactingMultiphaseParcelDevolatilisationModels.H"
#include "makeReactingMultiphaseParcelSurfaceReactionModels.H"

// MPPIC sub-models
#include "makeMPPICParcelDampingModels.H"
#include "makeMPPICParcelIsotropyModels.H"
#include "makeMPPICParcelPackingModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeReactingParcelCloudFunctionObjects(chargeCloud);

// Kinematic sub-models
makeThermoParcelForces(chargeCloud);
makeParcelDispersionModels(chargeCloud);
makeChargeParcelInjectionModels(chargeCloud);
makeParcelPatchInteractionModels(chargeCloud);
makeReactingMultiphaseParcelStochasticCollisionModels
(
    chargeCloud
);
makeReactingParcelSurfaceFilmModels(chargeCloud);

// Thermo sub-models
makeParcelHeatTransferModels(chargeCloud);

// Reacting sub-models
makeReactingMultiphaseParcelCompositionModels(chargeCloud);
makeReactingParcelPhaseChangeModels(chargeCloud);

// Reacting multiphase sub-models
makeReactingMultiphaseParcelDevolatilisationModels
(
    chargeCloud
);
makeReactingMultiphaseParcelSurfaceReactionModels
(
    chargeCloud
);

// MPPIC sub-models
makeMPPICParcelDampingModels(chargeCloud);
makeMPPICParcelIsotropyModels(chargeCloud);
makeMPPICParcelPackingModels(chargeCloud);

// ************************************************************************* //
