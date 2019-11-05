/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "vibrationShellModel.H"
#include "faMesh.H"
#include "fvMesh.H"
#include "fvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(vibrationShellModel, 0);

defineRunTimeSelectionTable(vibrationShellModel, dictionary);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


bool vibrationShellModel::read(const dictionary& dict)
{
    regionFaModel::read(dict);
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


vibrationShellModel::vibrationShellModel
(
    const word& modelType,
    const fvPatch& p,
    const dictionary& dict

)
:
    regionFaModel(p, "vibratingShell", modelType, dict, true),
    pName_(dict.get<word>("pa")),
    pa_(p.boundaryMesh().mesh().lookupObject<volScalarField>(pName_)),
    w_
    (
        IOobject
        (
            "ws_" + regionName_,
            p.boundaryMesh().mesh().time().timeName(),
            p.boundaryMesh().mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),
    a_
    (
        IOobject
        (
            "as_" + regionName_,
            p.boundaryMesh().mesh().time().timeName(),
            p.boundaryMesh().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimAcceleration, Zero)
    ),
    faOptions_(Foam::fa::options::New(p)),
    solid_(dict.subDict("solid"))
{
    if (!faOptions_.optionList::size())
    {
        Info << "No finite area options present" << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

vibrationShellModel::~vibrationShellModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void vibrationShellModel::preEvolveRegion()
{}

const Foam::volScalarField& vibrationShellModel::pa() const
{
    return pa_;
}


const Foam::areaScalarField& vibrationShellModel::w() const
{
    return w_;
}


const Foam::areaScalarField& vibrationShellModel::a() const
{
    return a_;
}


Foam::fa::options& vibrationShellModel::faOptions()
{
     return faOptions_;
}


const Foam::solidProperties& vibrationShellModel::solid() const
{
     return solid_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
