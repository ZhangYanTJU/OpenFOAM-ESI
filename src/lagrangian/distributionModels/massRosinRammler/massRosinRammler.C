/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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

#include "massRosinRammler.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace distributionModels
{
    defineTypeNameAndDebug(massRosinRammler, 0);
    addToRunTimeSelectionTable(distributionModel, massRosinRammler, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributionModels::massRosinRammler::massRosinRammler
(
    const dictionary& dict,
    Random& rndGen
)
:
    distributionModel(typeName, dict, rndGen),
    minValue_(distributionModelDict_.getScalar("minValue")),
    maxValue_(distributionModelDict_.getScalar("maxValue")),
    lambda_(distributionModelDict_.getCompat<scalar>("lambda", {{"d", 2012}})),
    n_(distributionModelDict_.getScalar("n"))
{
    if (lambda_ < VSMALL || n_ < VSMALL)
    {
        FatalErrorInFunction
            << "Scale/Shape parameter cannot be equal to or less than zero:"
            << "    lambda = " << lambda_
            << "    n = " << n_
            << exit(FatalError);
    }

    check();
}


Foam::distributionModels::massRosinRammler::massRosinRammler
(
    const massRosinRammler& p
)
:
    distributionModel(p),
    minValue_(p.minValue_),
    maxValue_(p.maxValue_),
    lambda_(p.lambda_),
    n_(p.n_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributionModels::massRosinRammler::sample() const
{
    /*scalar d;

    // Re-sample if the calculated d is out of the physical range
    do
    {
        // (YHD:Inverse of Eq. 10)
        const scalar a = scalar(3)/n_ + scalar(1);
        const scalar u = rndGen_.sample01<scalar>();
        const scalar x = invIncGamma(a, u);
        d = lambda_*pow(x, scalar(1)/n_);
    } while (d < minValue_ || d > maxValue_);

    return d;*/

    const scalar a = scalar(3)/n_ + scalar(1);
    const scalar cdfA = incGamma_P(a, pow(minValue_/lambda_, n_) );
    const scalar cdfB = incGamma_P(a, pow(maxValue_/lambda_, n_) );
    const scalar u = rndGen_.position<scalar>(cdfA, cdfB);
    const scalar x = invIncGamma(a, u);
    //return x;
    return lambda_*pow(x, scalar(1)/n_);
}


Foam::scalar Foam::distributionModels::massRosinRammler::minValue() const
{
    return minValue_;
}


Foam::scalar Foam::distributionModels::massRosinRammler::maxValue() const
{
    return maxValue_;
}


Foam::scalar Foam::distributionModels::massRosinRammler::meanValue() const
{
    // (YHD:Eqs. 11-12)
    return lambda_*tgamma(scalar(1)/n_ + scalar(1));
}


// ************************************************************************* //
