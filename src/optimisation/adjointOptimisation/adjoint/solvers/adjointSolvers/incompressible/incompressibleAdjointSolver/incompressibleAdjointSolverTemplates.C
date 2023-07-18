/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 PCOpt/NTUA
    Copyright (C) 2022-2023 FOSS GP
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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
tmp<surfaceInterpolationScheme<Type>>
Foam::incompressibleAdjointSolver::convectionScheme
(
    const word& varName,
    word phiName
) const
{
    const surfaceScalarField& phi = primalVars_.phi();
    if (phiName == word::null)
    {
        phiName = primalVars_.phiInst().name();
    }
    word divEntry("div(" + phiName + ',' + varName +')');
    ITstream& divScheme = mesh_.divScheme(divEntry);
    // Skip the first entry which might be 'bounded' or 'Gauss'.
    // If it is 'bounded', skip the second entry as well
    word discarded(divScheme);
    if (discarded == "bounded")
    {
        discarded = word(divScheme);
    }
    return surfaceInterpolationScheme<Type>::New(mesh_, phi, divScheme);
}


// ************************************************************************* //

