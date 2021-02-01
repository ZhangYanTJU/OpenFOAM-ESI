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

#include "normal.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace distributionModels
{
    defineTypeNameAndDebug(normal, 0);
    addToRunTimeSelectionTable(distributionModel, normal, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributionModels::normal::normal
(
    const dictionary& dict,
    Random& rndGen
)
:
    distributionModel(typeName, dict, rndGen),
    minValue_(distributionModelDict_.getScalar("minValue")),
    maxValue_(distributionModelDict_.getScalar("maxValue")),
    mu_(distributionModelDict_.getCompat<scalar>("mu",{{"expectation", 2012}})),
    sigmaSqr_
    (
        distributionModelDict_.getCompat<scalar>
        (
            "sigmaSqr",
            {{"variance", 2012}}
        )
    ),
    a_(0.147)
{
    if (mu_ < minValue_ || mu_ > maxValue_)
    {
        FatalErrorInFunction
            << "Expectation cannot be smaller than given minimum value, or "
            << "cannot be larger than given maximum value." << nl
            << "    mu = " << mu_ << nl
            << "    minValue = " << minValue_ << nl
            << "    maxValue = " << maxValue_
            << exit(FatalError);
    }

    check();
}


Foam::distributionModels::normal::normal(const normal& p)
:
    distributionModel(p),
    minValue_(p.minValue_),
    maxValue_(p.maxValue_),
    mu_(p.mu_),
    sigmaSqr_(p.sigmaSqr_),
    a_(p.a_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributionModels::normal::sample() const
{
    const scalar a = erf((minValue_ - mu_)/sigmaSqr_);
    const scalar b = erf((maxValue_ - mu_)/sigmaSqr_);

    const scalar u = rndGen_.sample01<scalar>();
    const scalar x = erfInv(u*(b - a) + a)*sigmaSqr_ + mu_;

    // Note: numerical approximation of the inverse function yields slight
    //       inaccuracies

    return min(max(x, minValue_), maxValue_);
}


Foam::scalar Foam::distributionModels::normal::minValue() const
{
    return minValue_;
}


Foam::scalar Foam::distributionModels::normal::maxValue() const
{
    return maxValue_;
}


Foam::scalar Foam::distributionModels::normal::meanValue() const
{
    return mu_;
}


Foam::scalar Foam::distributionModels::normal::erfInv(const scalar y) const
{
    scalar k = 2.0/(constant::mathematical::pi*a_) +  0.5*log(1.0 - y*y);
    scalar h = log(1.0 - y*y)/a_;
    scalar x = sqrt(-k + sqrt(k*k - h));
    if (y < 0.0)
    {
        x *= -1.0;
    }
    return x;
}


// ************************************************************************* //
