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

Class
    Foam::AMIFieldOpBase

Description
    Base class for AMI field operations

\*---------------------------------------------------------------------------*/

#ifndef AMIFieldOpsBase_H
#define AMIFieldOpsBase_H

#include "AMIInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class AMIFieldOpBase Declaration
\*---------------------------------------------------------------------------*/

class AMIFieldOpBase
{
protected:

    //- Reference to the AMI
    const AMIInterpolation& ami_;

    //- Interpolation to source; false=interpolate to target
    const bool toSource_;


public:

    //- Constructor
    AMIFieldOpBase(const AMIInterpolation& ami, const bool toSource)
    :
        ami_(ami),
        toSource_(toSource)
    {}


    // Member Functions

        //- Return the 'toSource' flag
        bool toSource() const
        {
            return toSource_;
        }

        //- Return the AMI addressing
        const labelListList& address() const
        {
            return toSource_ ? ami_.srcAddress() : ami_.tgtAddress();
        }

        //- Return the AMI weights
        const scalarListList& weights() const
        {
            return toSource_ ? ami_.srcWeights() : ami_.tgtWeights();
        }

        //- Return the AMI sum of weights
        const scalarList& weightsSum() const
        {
            return toSource_ ? ami_.srcWeightsSum() : ami_.tgtWeightsSum();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //