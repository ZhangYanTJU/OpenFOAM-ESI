/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2020,2025 OpenCFD Ltd.
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

{
    volScalarField& rDeltaT = trDeltaT.ref();

    const dictionary& pimpleDict = pimple.dict();

    // Maximum flow Courant number
    scalar maxCo(pimpleDict.get<scalar>("maxCo"));

    // Maximum time scale
    scalar maxDeltaT(pimpleDict.getOrDefault<scalar>("maxDeltaT", GREAT));

    // Smoothing parameter (0-1) when smoothing iterations > 0
    scalar rDeltaTSmoothingCoeff
    (
        pimpleDict.getOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.1)
    );

    // Damping coefficient (1-0)
    scalar rDeltaTDampingCoeff
    (
        pimpleDict.getOrDefault<scalar>("rDeltaTDampingCoeff", 1.0)
    );

    // Maximum change in cell temperature per iteration
    // (relative to previous value)
    scalar alphaTemp(pimpleDict.getOrDefault<scalar>("alphaTemp", 0.05));

    // Maximum change in cell concentration per iteration
    // (relative to reference value)
    scalar alphaY(pimpleDict.getOrDefault<scalar>("alphaY", 1.0));


    // The old reciprocal time scale field, with any damping factor
    tmp<volScalarField> rDeltaT0_damped;

    // Calculate damped value before applying any other changes
    if
    (
        rDeltaTDampingCoeff < 1
     && runTime.timeIndex() > runTime.startTimeIndex() + 1
    )
    {
        rDeltaT0_damped = (scalar(1) - rDeltaTDampingCoeff)*(rDeltaT);
    }


    Info<< "Time scales min/max:" << endl;

    // Flow time scale
    {
        rDeltaT.ref() =
        (
            fvc::surfaceSum(mag(phi))()()
           /((2*maxCo)*mesh.V()*rho())
        );

        // Limit the largest time scale (=> smallest reciprocal time)
        rDeltaT.clamp_min(1/maxDeltaT);

        auto limits = gMinMax(rDeltaT.primitiveField());
        limits.reset(1/(limits.max()+VSMALL), 1/(limits.min()+VSMALL));

        Info<< "    Flow        = "
            << limits.min() << ", " << limits.max() << endl;
    }

    // Heat release rate time scale
    if (alphaTemp < 1)
    {
        volScalarField::Internal rDeltaTT
        (
            mag(Qdot)/(alphaTemp*rho*thermo.Cp()*T)
        );

        rDeltaT.primitiveFieldRef().clamp_min(rDeltaTT);

        auto limits = gMinMax(rDeltaTT.field());
        limits.reset(1/(limits.max()+VSMALL), 1/(limits.min()+VSMALL));

        Info<< "    Temperature = "
            << limits.min() << ", " << limits.max() << endl;
    }

    // Reaction rate time scale
    if (alphaY < 1)
    {
        dictionary Yref(pimpleDict.subDict("Yref"));

        volScalarField::Internal rDeltaTY
        (
            IOobject
            (
                "rDeltaTY",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar(rDeltaT.dimensions(), Zero)
        );

        bool foundY = false;

        forAll(Y, i)
        {
            if (i != inertIndex && composition.active(i))
            {
                volScalarField& Yi = Y[i];

                if (Yref.found(Yi.name()))
                {
                    foundY = true;
                    const scalar Yrefi = Yref.get<scalar>(Yi.name());

                    rDeltaTY.field() = max
                    (
                        mag
                        (
                            reaction->R(Yi)().source()
                           /((Yrefi*alphaY)*(rho*mesh.V()))
                        ),
                        rDeltaTY
                    );
                }
            }
        }

        if (foundY)
        {
            rDeltaT.primitiveFieldRef().clamp_min(rDeltaTY);

            auto limits = gMinMax(rDeltaTY.field());
            limits.reset(1/(limits.max()+VSMALL), 1/(limits.min()+VSMALL));

            Info<< "    Composition = "
                << limits.min() << ", " << limits.max() << endl;
        }
        else
        {
            IOWarningIn(args.executable().c_str(), Yref)
                << "Cannot find any active species in Yref " << Yref
                << endl;
        }
    }

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    // Spatially smooth the time scale field
    if (rDeltaTSmoothingCoeff < 1)
    {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    }

    // Limit rate of change of time scale (=> smallest reciprocal time)
    // - reduce as much as required
    // - only increase at a fraction of old time scale
    if (rDeltaT0_damped)
    {
        rDeltaT.clamp_min(rDeltaT0_damped());
    }

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    auto limits = gMinMax(rDeltaT.field());
    limits.reset(1/(limits.max()+VSMALL), 1/(limits.min()+VSMALL));

    Info<< "    Overall     = "
        << limits.min() << ", " << limits.max() << endl;
}


// ************************************************************************* //
