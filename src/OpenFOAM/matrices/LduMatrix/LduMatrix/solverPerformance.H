/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2015 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

Typedef
    Foam::solverPerformance

Description
    SolverPerformance instantiated for a scalar

SourceFiles
    solverPerformance.C

\*---------------------------------------------------------------------------*/

#ifndef solverPerformance_H
#define solverPerformance_H

#include "SolverPerformance.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef SolverPerformance<scalar> solverPerformance;

    // Specialize max() member function for scalar
    template<>
    SolverPerformance<scalar> SolverPerformance<scalar>::max() const;

    // Specialize timing() member function for scalar
    template<>
    double SolverPerformance<scalar>::timing() const;

    // Specialize maxTiming() member function for scalar
    template<>
    double SolverPerformance<scalar>::maxTiming() const;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
