/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::pointPatchMapper

Description
    Mapping class for a pointPatchField.

SourceFiles
    pointPatchMapper.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_pointPatchMapper_H
#define Foam_pointPatchMapper_H

#include "pointMapper.H"
#include "pointPatchFieldMapper.H"
#include "pointPatch.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class pointPatch;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                       Class pointPatchMapper Declaration
\*---------------------------------------------------------------------------*/

class pointPatchMapper
:
    public pointPatchFieldMapper
{
    // Private Data

        //- Reference to patch
        const pointPatch& patch_;

        //- Reference to point mapper for internal field
        const morphFieldMapper& pointMapper_;

        //- Reference to mapping data
        const mapPolyMesh& mpm_;

        //- Size before mapping
        const label sizeBeforeMapping_;


    // Demand-driven Data

        mutable bool hasUnmapped_;

        //- Direct addressing (only one for of addressing is used)
        mutable std::unique_ptr<labelList> directAddrPtr_;

        //- Interpolated addressing (only one for of addressing is used)
        mutable std::unique_ptr<labelListList> interpAddrPtr_;

        //- Interpolation weights
        mutable std::unique_ptr<scalarListList> weightsPtr_;


    // Private Member Functions

        //- No copy construct
        pointPatchMapper(const pointPatchMapper&) = delete;

        //- No copy assignment
        void operator=(const pointPatchMapper&) = delete;


        //- Calculate addressing for mapping with inserted cells
        void calcAddressing() const;


public:

    // Static data members

    // Constructors

        //- Construct from mappers
        pointPatchMapper
        (
            const pointPatch& patch,
            const pointMapper& pointMap,
            const mapPolyMesh& mpm
        );


    //- Destructor
    virtual ~pointPatchMapper();


    // Member Functions

        //- Return size
        virtual label size() const
        {
            return patch_.size();
        }

        virtual bool hasUnmapped() const
        {
            return hasUnmapped_;
        }

        //- Return size of field before mapping
        virtual label sizeBeforeMapping() const
        {
            return sizeBeforeMapping_;
        }

        //- Is the mapping direct
        virtual bool direct() const
        {
            return pointMapper_.direct();
        }

        //- Return direct addressing
        virtual const labelUList& directAddressing() const;

        //- Return interpolated addressing
        virtual const labelListList& addressing() const;

        //- Return interpolaion weights
        virtual const scalarListList& weights() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
