/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2019 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "general.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace distributionModels
{
    defineTypeNameAndDebug(general, 0);
    addToRunTimeSelectionTable(distributionModel, general, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::distributionModels::general::initialise()
{
    static scalar eps = ROOTVSMALL;

    integral_.setSize(nEntries_);

    // Fill out the integral table (x, P(x<=0)) and calculate mean
    // For density function: P(x<=0) = int f(x) and mean = int x*f(x)
    // For cumulative function: mean = int 1-P(x<=0) = maxValue_ - int P(x<=0)
    integral_[0] = 0;

    for (label i = 1; i < nEntries_; ++i)
    {
        scalar k = (xy_[i][1] - xy_[i-1][1])/(xy_[i][0] - xy_[i-1][0] + eps);
        scalar d = xy_[i-1][1] - k*xy_[i-1][0];
        scalar y1 = xy_[i][0]*(0.5*k*xy_[i][0] + d);
        scalar y0 = xy_[i-1][0]*(0.5*k*xy_[i-1][0] + d);
        scalar area = y1 - y0;

        if (cumulative_)
        {
            integral_[i] = xy_[i][1];
            meanValue_ += area;
        }
        else
        {
            integral_[i] = area + integral_[i-1];

            y1 = sqr(xy_[i][0])*(1.0/3.0*k*xy_[i][0] + 0.5*d);
            y0 = sqr(xy_[i-1][0])*(1.0/3.0*k*xy_[i-1][0] + 0.5*d);
            meanValue_ += y1 - y0;
        }
    }

    // normalize the distribution
    scalar sumArea = integral_.last();

    for (label i = 0; i < nEntries_; ++i)
    {
        xy_[i][1] /= sumArea + eps;
        integral_[i] /= sumArea + eps;
    }

    meanValue_ /= sumArea;
    meanValue_ = cumulative_ ? (maxValue_ - meanValue_ + eps) : meanValue_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributionModels::general::general
(
    const dictionary& dict,
    Random& rndGen
)
:
    distributionModel(typeName, dict, rndGen),
    xy_(distributionModelDict_.lookup("distribution")),
    nEntries_(xy_.size()),
    minValue_(xy_[0][0]),
    maxValue_(xy_[nEntries_-1][0]),
    meanValue_(0),
    integral_(nEntries_),
    cumulative_(distributionModelDict_.getOrDefault("cumulative", false))
{
    check();

    // Additional sanity checks
    if (cumulative_ && xy_[0][1] != 0)
    {
        FatalErrorInFunction
            << type() << "distribution: "
            << "Elements in the second column for cumulative "
            << "distribution functions must start from zero." << nl
            << "First element = " << xy_[0][1]
            << exit(FatalError);
    }

    for (label i = 0; i < nEntries_; ++i)
    {
        if (i > 0 && xy_[i][0] <= xy_[i-1][0])
        {
            FatalErrorInFunction
                << type() << "distribution: "
                << "Elements in the first column must "
                << "be specified in an ascending order." << nl
                << "Please see the row i = " << i << nl
                << "xy[i-1] = " << xy_[i-1] << nl
                << "xy[i] = " << xy_[i]
                << exit(FatalError);
        }

        if (xy_[i][0] < 0 || xy_[i][1] < 0)
        {
            FatalErrorInFunction
                << type() << "distribution: "
                << "Input pairs cannot contain any negative element." << nl
                << "Please see the row i = " << i << nl
                << "xy[i] = " << xy_[i]
                << exit(FatalError);
        }

        if (cumulative_ && i > 0 && xy_[i][1] < xy_[i-1][1])
        {
            FatalErrorInFunction
                << type() << "distribution: "
                << "Elements in the second column for cumulative "
                << "distribution functions must be non-decreasing." << nl
                << "Please see the row i = " << i << nl
                << "xy[i-1] = " << xy_[i-1] << nl
                << "xy[i] = " << xy_[i]
                << exit(FatalError);
        }
    }

    initialise();
}


Foam::distributionModels::general::general
(
    const UList<scalar>& sampleData,
    const scalar binWidth,
    Random& rndGen
)
:
    distributionModel(typeName, dictionary::null, rndGen),
    xy_(),
    nEntries_(0),
    minValue_(0),
    maxValue_(0),
    meanValue_(0),
    integral_(),
    cumulative_(false)
{
    scalar minValue = GREAT;
    scalar maxValue = -GREAT;
    forAll(sampleData, i)
    {
        minValue = min(minValue, sampleData[i]);
        maxValue = max(maxValue, sampleData[i]);
    }

    label bin0 = floor(minValue/binWidth);
    label bin1 = ceil(maxValue/binWidth);
    nEntries_ = bin1 - bin0;

    if (nEntries_ == 0)
    {
        WarningInFunction
            << "Data cannot be binned - zero bins generated" << nl
            << "   Bin width   : " << binWidth << nl
            << "   Sample data : " << sampleData
            << endl;

        return;
    }

    xy_.setSize(nEntries_);

    // Populate bin boundaries and initialise occurrences
    for (label bini = 0; bini < nEntries_; ++bini)
    {
        xy_[bini][0] = (bin0 + bini)*binWidth;
        xy_[bini][1] = 0;
    }

    // Populate occurrences
    forAll(sampleData, i)
    {
        label bini = floor(sampleData[i]/binWidth) - bin0;
        xy_[bini][1]++;
    }

    initialise();
}


Foam::distributionModels::general::general(const general& p)
:
    distributionModel(p),
    xy_(p.xy_),
    nEntries_(p.nEntries_),
    minValue_(p.minValue_),
    maxValue_(p.maxValue_),
    meanValue_(p.meanValue_),
    integral_(p.integral_),
    cumulative_(p.cumulative_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributionModels::general::sample() const
{
    scalar y = rndGen_.sample01<scalar>();

    // Find the interval where y is in the table
    label n = 1;
    while (integral_[n] <= y)
    {
        n++;
    }

    scalar k = (xy_[n][1] - xy_[n-1][1])/(xy_[n][0] - xy_[n-1][0]);
    scalar d = xy_[n-1][1] - k*xy_[n-1][0];

    if (cumulative_)
    {
        return (y - d)/k;
    }

    scalar alpha = y + xy_[n-1][0]*(0.5*k*xy_[n-1][0] + d) - integral_[n-1];
    scalar x = 0;

    // If k is small it is a linear equation, otherwise it is of second order
    if (mag(k) > SMALL)
    {
        scalar p = 2.0*d/k;
        scalar q = -2.0*alpha/k;
        scalar sqrtEr = sqrt(0.25*p*p - q);

        scalar x1 = -0.5*p + sqrtEr;
        scalar x2 = -0.5*p - sqrtEr;
        if ((x1 >= xy_[n-1][0]) && (x1 <= xy_[n][0]))
        {
            x = x1;
        }
        else
        {
            x = x2;
        }
    }
    else
    {
        x = alpha/d;
    }

    return x;
}


Foam::scalar Foam::distributionModels::general::minValue() const
{
    return minValue_;
}


Foam::scalar Foam::distributionModels::general::maxValue() const
{
    return maxValue_;
}


Foam::scalar Foam::distributionModels::general::meanValue() const
{
    return meanValue_;
}


void Foam::distributionModels::general::readData(Istream& is)
{
//    distributionModel::readData(is);
    is  >> xy_;
    initialise();
}


void Foam::distributionModels::general::writeData(Ostream& os) const
{
//    distributionModel::writeData(os);
    os  << xy_;
}


Foam::dictionary Foam::distributionModels::general::writeDict
(
    const word& dictName
) const
{
    // dictionary dict = distributionModel::writeDict(dictName);
    dictionary dict(dictName);
    dict.add("x", x());
    dict.add("y", y());

    return dict;
}


void Foam::distributionModels::general::readDict(const dictionary& dict)
{
    // distributionModel::readDict(dict);
    List<scalar> x(dict.lookup("x"));
    List<scalar> y(dict.lookup("y"));

    xy_.setSize(x.size());
    forAll(xy_, i)
    {
        xy_[i][0] = x[i];
        xy_[i][1] = y[i];
    }

    initialise();
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::distributionModels::general::x() const
{
    auto tx = tmp<scalarField>::New(xy_.size());
    auto& xi = tx.ref();
    forAll(xy_, i)
    {
        xi[i] = xy_[i][0];
    }

    return tx;
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::distributionModels::general::y() const
{
    auto ty = tmp<scalarField>::New(xy_.size());
    auto& yi = ty.ref();
    forAll(xy_, i)
    {
        yi[i] = xy_[i][1];
    }

    return ty;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const distributionModels::general& b
)
{
    os.check(FUNCTION_NAME);

    b.writeData(os);
    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, distributionModels::general& b)
{
    is.check(FUNCTION_NAME);

    b.readData(is);
    return is;
}


// ************************************************************************* //
