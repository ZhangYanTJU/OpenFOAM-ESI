/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2024 OpenCFD Ltd.
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

#ifndef makeThermoParcelCloudFunctionObjects_H
#define makeThermoParcelCloudFunctionObjects_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FaceInteraction.H"
#include "FacePostProcessing.H"
#include "FreezeParticles.H"
#include "ParticleCollector.H"
#include "ParticleErosion.H"
#include "ParticleTracks.H"
#include "ParticleTrap.H"
#include "ParticleZoneInfo.H"
#include "ParticleHistogram.H"
#include "ParticlePostProcessing.H"
#include "PatchCollisionDensity.H"
#include "PatchInteractionFields.H"
#include "RemoveParcels.H"
#include "VoidFraction.H"
#include "NusseltNumber.H"
#include "HeatTransferCoeff.H"
#include "ThermoReynoldsNumber.H"
#include "ParticleDose.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeThermoParcelCloudFunctionObjects(CloudType)                        \
                                                                               \
    makeCloudFunctionObject(CloudType);                                        \
                                                                               \
    makeCloudFunctionObjectType(FaceInteraction, CloudType);                   \
    makeCloudFunctionObjectType(FacePostProcessing, CloudType);                \
    makeCloudFunctionObjectType(FreezeParticles, CloudType);                   \
    makeCloudFunctionObjectType(ParticleCollector, CloudType);                 \
    makeCloudFunctionObjectType(ParticleErosion, CloudType);                   \
    makeCloudFunctionObjectType(ParticleTracks, CloudType);                    \
    makeCloudFunctionObjectType(ParticleTrap, CloudType);                      \
    makeCloudFunctionObjectType(ParticleZoneInfo, CloudType);                  \
    makeCloudFunctionObjectType(ParticleHistogram, CloudType);                 \
    makeCloudFunctionObjectType(ParticlePostProcessing, CloudType);            \
    makeCloudFunctionObjectType(PatchCollisionDensity, CloudType);             \
    makeCloudFunctionObjectType(PatchInteractionFields, CloudType);            \
    makeCloudFunctionObjectType(RemoveParcels, CloudType);                     \
    makeCloudFunctionObjectType(VoidFraction, CloudType);                      \
    makeCloudFunctionObjectType(NusseltNumber, CloudType);                     \
    makeCloudFunctionObjectType(HeatTransferCoeff, CloudType);                 \
    makeCloudFunctionObjectType(ThermoReynoldsNumber, CloudType);              \
    makeCloudFunctionObjectType(ParticleDose, CloudType);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
