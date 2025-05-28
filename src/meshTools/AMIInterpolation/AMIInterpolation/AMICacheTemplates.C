/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenCFD Ltd.
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class EvalFunction>
bool Foam::AMICache::apply(List<Type>& result, const EvalFunction& eval) const
{
    if (applyLower())
    {
        eval(result, index0_);
        return true;
    }
    else if (applyUpper())
    {
        eval(result, index1_);
        return true;
    }
    else if (applyInterpolate())
    {
        List<Type> r0(result);
        eval(r0, index0_);
        List<Type> r1(result);
        eval(r1, index1_);

        //result = (r1 - r0)*interpWeight_ + r0;
        forAll(result, i)
        {
            result[i] = lerp(r0[i], r1[i], interpWeight_);
        }
        return true;
    }

    return false;
}


// ************************************************************************* //
