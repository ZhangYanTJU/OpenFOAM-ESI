/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "tabulated6DoFMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "Tuple2.H"
#include "IFstream.H"
#include "interpolateSplineXY.H"
#include "interpolateXY.H"
#include "unitConversion.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(tabulated6DoFMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        tabulated6DoFMotion,
        dictionary
    );
}
}


const Foam::Enum
<
    Foam::solidBodyMotionFunctions::tabulated6DoFMotion::interpolationType
>
Foam::solidBodyMotionFunctions::tabulated6DoFMotion::interpolationTypeNames
({
    { interpolationType::SPLINE, "spline" },
    { interpolationType::LINEAR, "linear" }
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::tabulated6DoFMotion::tabulated6DoFMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime)
{
    read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::tabulated6DoFMotion::transformation() const
{
    scalar t = time_.value();

    if (t < times_[0])
    {
        FatalErrorInFunction
            << "current time (" << t
            << ") is less than the minimum in the data table ("
            << times_[0] << ')'
            << exit(FatalError);
    }

    if (t > times_.last())
    {
        FatalErrorInFunction
            << "current time (" << t
            << ") is greater than the maximum in the data table ("
            << times_.last() << ')'
            << exit(FatalError);
    }

    translationRotationVectors TRV(Zero, Zero);
    switch (interpolator_)
    {
        case interpolationType::SPLINE:
        {
            TRV = interpolateSplineXY(t, times_, values_);
            break;
        }
        case interpolationType::LINEAR:
        {
            TRV = interpolateXY(t, times_, values_);
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unrecognised 'interpolationScheme' option: "
                << interpolationTypeNames[interpolator_]
                << exit(FatalError);
            break;
        }
    }

    // Convert the rotational motion from deg to rad
    TRV[1] *= degToRad();

    quaternion R(quaternion::XYZ, TRV[1]);
    septernion TR(septernion(-CofG_ + -TRV[0])*R*septernion(CofG_));

    DebugInFunction << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::tabulated6DoFMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    // If the timeDataFileName has changed read the file

    fileName newTimeDataFileName
    (
        SBMFCoeffs_.get<fileName>("timeDataFileName").expand()
    );

    if (newTimeDataFileName != timeDataFileName_)
    {
        timeDataFileName_ = newTimeDataFileName;

        IFstream dataStream(timeDataFileName_);

        if (dataStream.good())
        {
            List<Tuple2<scalar, translationRotationVectors>> timeValues
            (
                dataStream
            );

            times_.setSize(timeValues.size());
            values_.setSize(timeValues.size());

            forAll(timeValues, i)
            {
                times_[i] = timeValues[i].first();
                values_[i] = timeValues[i].second();
            }
        }
        else
        {
            FatalErrorInFunction
                << "Cannot open time data file " << timeDataFileName_
                << exit(FatalError);
        }
    }

    SBMFCoeffs_.readEntry("CofG", CofG_);

    interpolator_ =
        interpolationTypeNames.getOrDefault
        (
            "interpolationScheme",
            SBMFCoeffs_,
            interpolationType::SPLINE
        );

    return true;
}


// ************************************************************************* //
