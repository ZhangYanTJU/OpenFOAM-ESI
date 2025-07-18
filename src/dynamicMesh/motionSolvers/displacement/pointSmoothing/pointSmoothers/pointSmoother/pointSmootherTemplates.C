/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024-2025 OpenCFD Ltd.
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

#include "pointSmoother.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template <class weightType>
void Foam::pointSmoother::reset
(
    const labelList& facesToMove,
    Field<weightType>& weights,
    vectorField& pointDisplacement,
    const bool resetInternalFaces
) const
{
    bitSet resetPoints
    (
        pointsToMove(facesToMove, resetInternalFaces)
    );

    for (const label pointi: resetPoints)
    {
        weights[pointi] = pTraits<weightType>::zero;
        pointDisplacement[pointi] = vector::zero;
    }
}


template <class weightType>
void Foam::pointSmoother::average
(
    const labelList& facesToMove,
    Field<weightType>& weights,
    vectorField& pointDisplacement
) const
{
    syncTools::syncPointList
    (
        mesh(),
        weights,
        plusEqOp<weightType>(),
        pTraits<weightType>::zero
    );

    syncTools::syncPointList
    (
        mesh(),
        pointDisplacement,
        plusEqOp<vector>(),
        vector::zero
    );

    bitSet averagePoints
    (
        pointsToMove(facesToMove, true)
    );

    for (const label pointi : averagePoints)
    {
        if (weights[pointi] != pTraits<weightType>::zero)
        {
            pointDisplacement[pointi] /= weights[pointi];
        }
    }
}


// ************************************************************************* //
