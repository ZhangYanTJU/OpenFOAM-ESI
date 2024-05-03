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
    Foam::LowWeightCorrection

Description
    Wrapper around the AMI field op

\*---------------------------------------------------------------------------*/

#ifndef LowWeightCorrection_H
#define LowWeightCorrection_H

#include "lowWeightCorrectionBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class LowWeightCorrection Declaration
\*---------------------------------------------------------------------------*/

template<class AMIOpType, class CombineOp>
class LowWeightCorrection
:
    public lowWeightCorrectionBase,
    public AMIOpType
{
    typedef typename AMIOpType::value_type value_type;

public:

    //- Constructor
    LowWeightCorrection
    (
        const AMIInterpolation& ami,
        const bool toSource,
        const lowWeightCorrectionBase::option& opt =
            lowWeightCorrectionBase::option::NONE
    )
    :
        lowWeightCorrectionBase(opt),
        AMIOpType(ami, toSource)
    {}


    // Member Functions

        //- Helper function to ensure default value field is valid
        bool validDefaults(const UList<value_type>& defaultValues) const
        {
            const auto& ami = AMIOpType::ami_;

            if (ami.lowWeightCorrection() > 0)
            {
                const auto& weightsSum = this->weightsSum();

                if (defaultValues.size() != weightsSum.size())
                {
                    FatalErrorInFunction
                        << "Employing default values when sum of weights "
                        << "falls below " << ami.lowWeightCorrection()
                        << " but number of default values is not equal to "
                        << "addressing size" << nl
                        << "    default values size = "
                        << defaultValues.size() << nl
                        << "    addressing size     = "
                        << weightsSum.size() << nl
                        << abort(FatalError);
                }

                return true;
            }

            return false;
        }


        //- Correction operator
        void operator()
        (
            List<value_type>& result,
            const UList<value_type>& fld,
            const UList<value_type>& defaultValues
        ) const
        {
            // Apply AMI interpolator
            AMIOpType::operator()(result, fld, UList<value_type>::null());

            if (!validDefaults(defaultValues)) return;

            switch (opt_)
            {
                case option::NONE:
                {
                    break;
                }
                case option::ASSIGN:
                {
                    const auto& ami = AMIOpType::ami_;
                    const auto& weightsSum = this->weightsSum();

                    forAll(result, facei)
                    {
                        if (weightsSum[facei] < ami.lowWeightCorrection())
                        {
                            // Assign default value
                            result[facei] = defaultValues[facei];
                        }
                    }

                    break;
                }
                case option::BLEND:
                {
                    const auto& ami = AMIOpType::ami_;
                    const auto& weightsSum = this->weightsSum();

                    forAll(result, facei)
                    {
                        if (weightsSum[facei] < ami.lowWeightCorrection())
                        {
                            // Already have 'weightSum' contribution in result
                            // - blend 'missing' contribution into result
                            const scalar w = min(1, max(0, weightsSum[facei]));
                            result[facei] += (1 - w)*defaultValues[facei];
                        }
                    }

                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Unhandled enumeration " << static_cast<int>(opt_)
                        << abort(FatalError);
                }
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //