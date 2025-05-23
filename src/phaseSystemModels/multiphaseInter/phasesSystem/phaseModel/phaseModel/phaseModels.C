/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "makePhaseTypes.H"

#include "PurePhaseModel.H"
#include "MultiComponentPhaseModel.H"
#include "MovingPhaseModel.H"
#include "StaticPhaseModel.H"

#include "rhoThermo.H"
#include "solidThermo.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePhaseTypes
(
    MovingPhaseModel,
    PurePhaseModel,
    multiphaseInter::phaseModel,
    rhoThermo,
    pureMovingPhaseModel // Name of the phase type
);

makePhaseTypes
(
    StaticPhaseModel,
    PurePhaseModel,
    multiphaseInter::phaseModel,
    rhoThermo,
    pureStaticPhaseModel
);

makePhaseTypes
(
    StaticPhaseModel,
    PurePhaseModel,
    multiphaseInter::phaseModel,
    solidThermo,
    pureStaticSolidPhaseModel
);

makePhaseTypes
(
    MovingPhaseModel,
    PurePhaseModel,
    multiphaseInter::phaseModel,
    solidThermo,
    pureMovingSolidPhaseModel
);

makePhaseTypes
(
    MovingPhaseModel,
    MultiComponentPhaseModel,
    multiphaseInter::phaseModel,
    rhoReactionThermo,
    multiComponentMovingPhaseModel
);


// ************************************************************************* //
