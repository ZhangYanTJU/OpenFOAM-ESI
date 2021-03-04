/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "KTS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace XiGModels
{
    defineTypeNameAndDebug(KTS, 0);
    addToRunTimeSelectionTable(XiGModel, KTS, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiGModels::KTS::KTS
(
    const dictionary& XiGProperties,
    const word& modelType,
    const psiuReactionThermo& thermo,
    const compressible::RASModel& turbulence,
    const volScalarField& Su
)
:
    XiGModel(XiGProperties, modelType, thermo, turbulence, Su),
    GEtaCoef_(XiGModelCoeffs_.get<scalar>("GEtaCoef"))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::XiGModels::KTS::G() const
{
    volScalarField up(sqrt((2.0/3.0)*turbulence_.k()));
    const volScalarField& epsilon = turbulence_.epsilon();

    volScalarField tauEta(sqrt(mag(thermo_.muu()/(thermo_.rhou()*epsilon))));

    return (GEtaCoef_/tauEta);
}


Foam::tmp<Foam::volScalarField> Foam::XiGModels::KTS::Db() const
{
    const objectRegistry& db = Su_.db();
    const volScalarField& Db1 = db.lookupObject<volScalarField>("Db");
    //return turbulence_.muEff();
    return Db1;
}


bool Foam::XiGModels::KTS::read(const dictionary& XiGProperties)
{
    XiGModel::read(XiGProperties);

    XiGModelCoeffs_.readEntry("GEtaCoef", GEtaCoef_);

    return true;
}


// ************************************************************************* //
