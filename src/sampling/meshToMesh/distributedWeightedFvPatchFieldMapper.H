/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
    Copyright (C) 2015-2019 OpenCFD Ltd.
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
    Foam::distributedWeightedFvPatchFieldMapper

Description
    FieldMapper with weighted mapping from (optionally remote) quantities.

\*---------------------------------------------------------------------------*/

#ifndef distributedWeightedFvPatchFieldMapper_H
#define distributedWeightedFvPatchFieldMapper_H

#include "fvPatchFieldMapper.H"
#include "mapDistributeBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class distributedWeightedFvPatchFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class distributedWeightedFvPatchFieldMapper
:
    public fvPatchFieldMapper
{
    const mapDistributeBase* distMapPtr_;

    const labelListList& addressing_;

    const scalarListList& weights_;

    bool hasUnmapped_;

public:

    // Constructors

        //- Construct given addressing
        distributedWeightedFvPatchFieldMapper
        (
            const label singlePatchProc,            // global consistent
            const label comm,                       // only on processors with
                                                    // patch faces
            const mapDistributeBase* distMapPtr,    //       ,,
            const labelListList& addressing,
            const scalarListList& weights
        )
        :
            distMapPtr_(distMapPtr),
            addressing_(addressing),
            weights_(weights),
            hasUnmapped_(false)
        {
            for (const labelList& addr : addressing)
            {
                if (addr.empty())
                {
                    hasUnmapped_ = true;
                    break;
                }
            }


            // singlePatchProc : >= 0
            //  - comm not relevant
            //  - distMap should not be set
            // singlePatchProc : -1
            //  - if (comm != -1) we need a distMap
            bool error = false;
            if (singlePatchProc >= 0)
            {
                error = distMapPtr_;
            }
            else
            {
                // Need map (unless I don't have communicator)
                if (comm >= 0)
                {
                    error = !distMapPtr_;
                }
                else
                {
                    error = distMapPtr_;
                }
            }
            if (error)
            {
                FatalErrorInFunction
                    << "Supply a mapDistributeBase if and only if "
                    << "singlePatchProc is -1"
                    << " singlePatchProc:" << singlePatchProc
                    << " distMapPtr_:" << (distMapPtr_ != nullptr)
                    << exit(FatalError);
            }
        }

    //- Destructor
    virtual ~distributedWeightedFvPatchFieldMapper() = default;


    // Member Functions

        virtual label size() const
        {
            if (distributed())
            {
                return distributeMap().constructSize();
            }
            else
            {
                return addressing().size();
            }
        }

        virtual bool direct() const
        {
            return false;
        }

        virtual bool distributed() const
        {
            return distMapPtr_;
        }

        virtual const mapDistributeBase& distributeMap() const
        {
            if (!distMapPtr_)
            {
                FatalErrorInFunction
                    << "Cannot ask for distributeMap on a non-distributed"
                    << " mapper" << exit(FatalError);
            }
            return *distMapPtr_;
        }

        virtual bool hasUnmapped() const
        {
            return hasUnmapped_;
        }

        virtual const labelListList& addressing() const
        {
            return addressing_;
        }

        virtual const scalarListList& weights() const
        {
            return weights_;
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
