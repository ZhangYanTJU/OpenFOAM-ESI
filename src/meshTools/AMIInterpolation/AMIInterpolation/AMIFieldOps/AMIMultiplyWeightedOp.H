/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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

#ifndef AMIMultiplyWeightedOp_H
#define AMIMultiplyWeightedOp_H

#include "AMIFieldOpBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class AMIMultiplyWeightedOp Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class CombineOp = plusEqOp<Type>>
class AMIMultiplyWeightedOp
:
    public AMIFieldOpBase
{
    // Private Data

        //- Combine operator, e.g. plusEqOp<Type>()
        //  Note: must be null-constructed
        CombineOp cop_;


public:

    typedef Type value_type;


    //- Constructor
    AMIMultiplyWeightedOp
    (
        const AMIInterpolation& ami,
        const bool toSource
    )
    :
        AMIFieldOpBase(ami, toSource),
        cop_()
    {}


    //- Multiply-weighted op
    //  Note: does not use default values
    void operator()
    (
        List<value_type>& result,
        const UList<value_type>& fld,
        const UList<value_type>& /* unused defaultValues */
    ) const
    {
        const auto& address = this->address();
        const auto& weights = this->weights();

        forAll(result, facei)
        {
            const labelList& faceSlots = address[facei];
            const scalarList& faceWeights = weights[facei];

            forAll(faceSlots, i)
            {
                if (i == 0)
                {
                    result[facei] = faceWeights[i]*fld[faceSlots[i]];
                }
                else
                {
                    cop_(result[facei], faceWeights[i]*fld[faceSlots[i]]);
                }
            }
        }
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //