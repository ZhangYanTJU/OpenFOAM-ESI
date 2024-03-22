/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "evapotranspirationHeatTransferModel.H"
#include "radiationModel.H"
#include "solarLoadBase.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(evapotranspirationHeatTransferModel, 0);
    defineRunTimeSelectionTable
    (
        evapotranspirationHeatTransferModel,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::volScalarField& Foam::evapotranspirationHeatTransferModel::getOrReadField
(
    const word& fieldName
) const
{
    auto* ptr = mesh_.getObjectPtr<volScalarField>(fieldName);

    if (!ptr)
    {
        ptr = new volScalarField
        (
            IOobject
            (
                fieldName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        );
        mesh_.objectRegistry::store(ptr);
    }

    return *ptr;
}


Foam::scalar Foam::evapotranspirationHeatTransferModel::q() const
{
    // Retrieve solar-load information
    const auto& base =
        mesh().lookupObject<radiation::solarLoadBase>("solarLoadBase");
    const solarCalculator& solar = base.solarCalculatorRef();

    // Retrieve internal & boundary faces directly hit by solar rays
    const faceShading& shade = base.faceShadingRef();
    const labelList& hitFacesId = shade.rayStartFaces();

    // Retrieve face zone information
    const faceZone& zone = mesh_.faceZones()[faceZoneId_];
    const vectorField& faceNormals = zone().faceNormals();
    const scalarField& faceAreas = zone().magFaceAreas();

    // Retrieve direct solar load [W/m^2]
    const vector directSolarRad(solar.directSolarRad()*solar.direction());

    // Identify zone faces within mesh
    bitSet isZoneFace(mesh_.nFaces());
    isZoneFace.set(zone);

    // Calculate area-averaged incident solar load
    // Assuming hit faces are updated by solarLoad
    scalar q = 0;
    scalar totalFaceArea = 0;

    for (const label facei : hitFacesId)
    {
        if (isZoneFace[facei])
        {
            const label localFacei = zone.whichFace(facei);

            const vector& faceNormal = faceNormals[localFacei];
            const scalar faceArea = faceAreas[localFacei];

            const scalar qIncident = directSolarRad & faceNormal;

            q += qIncident*faceArea;
            totalFaceArea += faceArea;
        }
    }
    reduce(q, sumOp<scalar>());
    reduce(totalFaceArea, sumOp<scalar>());
    q /= totalFaceArea;


    // Sum diffusive solar loads
    switch(solar.sunLoadModel())
    {
        case solarCalculator::mSunLoadFairWeatherConditions:
        case solarCalculator::mSunLoadTheoreticalMaximum:
        {
            // Calculate diffusive radiance
            // contribution from sky and ground
            tmp<scalarField> tdiffuseSolarRad =
                solar.diffuseSolarRad(faceNormals);
            const scalarField& diffuseSolarRad = tdiffuseSolarRad.cref();

            // Calculate area-averaged diffusive solar load
            scalar meanDiffuseSolarRad = 0;
            forAll(faceAreas, i)
            {
                meanDiffuseSolarRad += diffuseSolarRad[i]*faceAreas[i];
            }
            meanDiffuseSolarRad /= totalFaceArea;

            q += meanDiffuseSolarRad;
        }
        break;

        case solarCalculator::mSunLoadConstant:
        case solarCalculator::mSunLoadTimeDependent:
        {
            q += solar.diffuseSolarRad();
        }
        break;
    }

    return q;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::evapotranspirationHeatTransferModel::evapotranspirationHeatTransferModel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::evapotranspirationHeatTransferModel::read(const dictionary& dict)
{
    faceZoneId_ = mesh_.faceZones().findZoneID(dict.get<word>("faceZone"));

    if (faceZoneId_ < 0)
    {
        const word faceZoneName(mesh_.faceZones()[faceZoneId_].name());

        FatalIOErrorInFunction(dict)
            << type() << ' ' << typeName << ": "
            << "    No matching face zone: " << faceZoneName  << nl
            << "    Known face zones: "
            << flatOutput(mesh_.faceZones().names()) << nl
            << exit(FatalIOError);

        return false;
    }

    return true;
}


// ************************************************************************* //
