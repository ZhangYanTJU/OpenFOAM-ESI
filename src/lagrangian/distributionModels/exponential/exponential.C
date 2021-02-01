/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "exponential.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace distributionModels
{
    defineTypeNameAndDebug(exponential, 0);
    addToRunTimeSelectionTable(distributionModel, exponential, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributionModels::exponential::exponential
(
    const dictionary& dict,
    Random& rndGen
)
:
    distributionModel(typeName, dict, rndGen),
    minValue_(distributionModelDict_.getScalar("minValue")),
    maxValue_(distributionModelDict_.getScalar("maxValue")),
    lambda_(distributionModelDict_.getScalar("lambda"))
{
    if (lambda_ < VSMALL)
    {
        FatalErrorInFunction
            << "Rate parameter cannot be equal to or less than zero:"
            << "    lambda = " << lambda_
            << exit(FatalError);
    }

    check();
}


Foam::distributionModels::exponential::exponential(const exponential& p)
:
    distributionModel(p),
    minValue_(p.minValue_),
    maxValue_(p.maxValue_),
    lambda_(p.lambda_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributionModels::exponential::sample() const
{
    const scalar u = rndGen_.sample01<scalar>();
    const scalar qMin = exp(-lambda_*minValue_);
    const scalar qMax = exp(-lambda_*maxValue_);
    return -(scalar(1)/lambda_)*log(qMin + u*(qMax - qMin));
}


Foam::scalar Foam::distributionModels::exponential::minValue() const
{
    return minValue_;
}


Foam::scalar Foam::distributionModels::exponential::maxValue() const
{
    return maxValue_;
}


Foam::scalar Foam::distributionModels::exponential::meanValue() const
{
    return scalar(1)/lambda_;
}


// ************************************************************************* //
