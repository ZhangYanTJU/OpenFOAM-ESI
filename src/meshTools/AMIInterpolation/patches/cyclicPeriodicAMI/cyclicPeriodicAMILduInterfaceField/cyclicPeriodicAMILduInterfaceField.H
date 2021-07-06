/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::cyclicPeriodicAMILduInterfaceField

Description
    Abstract base class for cyclic periodic AMI coupled interfaces

SourceFiles
    cyclicPeriodicAMILduInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicPeriodicAMILduInterfaceField_H
#define cyclicPeriodicAMILduInterfaceField_H

#include "cyclicAMILduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class cyclicPeriodicAMILduInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class cyclicPeriodicAMILduInterfaceField
:
    public cyclicAMILduInterfaceField
{

public:

    //- Runtime type information
    TypeName("cyclicPeriodicAMILduInterfaceField");


    // Constructors

        //- Construct null
        cyclicPeriodicAMILduInterfaceField()
        {}


    //- Destructor
    virtual ~cyclicPeriodicAMILduInterfaceField();


    // Member Functions

        //- Transform given patch field
        template<class Type>
        void transformCoupleField(Field<Type>& f) const;

        //- Transform given patch internal field
        void transformCoupleField
        (
            solveScalarField& psiInternal,
            const direction cmpt
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tensorField.H"

template<class Type>
void Foam::cyclicPeriodicAMILduInterfaceField::transformCoupleField
(
    Field<Type>& f
) const
{
    cyclicAMILduInterfaceField::transformCoupleField(f);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
