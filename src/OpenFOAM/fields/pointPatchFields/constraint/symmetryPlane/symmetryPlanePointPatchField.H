/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
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
    Foam::symmetryPlanePointPatchField

Description
    A symmetry-plane boundary condition for pointField

SourceFiles
    symmetryPlanePointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef symmetryPlanePointPatchField_H
#define symmetryPlanePointPatchField_H

#include "basicSymmetryPointPatchField.H"
#include "symmetryPlanePointPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class symmetryPlanePointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class symmetryPlanePointPatchField
:
    public basicSymmetryPointPatchField<Type>
{
    // Private data

        //- Local reference cast into the symmetryPlane patch
        const symmetryPlanePointPatch& symmetryPlanePatch_;


public:

    //- Runtime type information
    TypeName(symmetryPlanePointPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        symmetryPlanePointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        symmetryPlanePointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch field onto a new patch
        symmetryPlanePointPatchField
        (
            const symmetryPlanePointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        symmetryPlanePointPatchField
        (
            const symmetryPlanePointPatchField<Type>&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Return a clone
        virtual autoPtr<pointPatchField<Type>> clone() const
        {
            return pointPatchField<Type>::Clone(*this);
        }

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<Type>> clone
        (
            const DimensionedField<Type, pointMesh>& iF
        ) const
        {
            return pointPatchField<Type>::Clone(*this, iF);
        }


    // Member functions

        //- Return the constraint type this pointPatchField implements
        virtual const word& constraintType() const
        {
            return symmetryPlanePointPatch::typeName;
        }

        //- Update the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::commsTypes::buffered
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "symmetryPlanePointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
