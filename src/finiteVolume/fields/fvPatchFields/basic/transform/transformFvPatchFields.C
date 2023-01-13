/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "transformFvPatchFields.H"
#include "fvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makePatchFieldTypeNames(transform);
    makePatchFieldTypeName(label, transform);
}


// * * * * * * * * * * * * * * * Specialisations * * * * * * * * * * * * * * //

template<>
Foam::tmp<Foam::Field<Foam::label>>
Foam::transformFvPatchField<Foam::label>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    // TBD: Treat like zero-gradient
    return tmp<Field<label>>::New(this->size(), label(1));
}


template<>
Foam::tmp<Foam::Field<Foam::label>>
Foam::transformFvPatchField<Foam::label>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    // TBD: Treat like zero-gradient
    return tmp<Field<label>>::New(this->size(), Zero);
}


template<>
Foam::tmp<Foam::Field<Foam::label>>
Foam::transformFvPatchField<Foam::label>::gradientInternalCoeffs() const
{
    // TBD: Treat like zero-gradient
    return tmp<Field<label>>::New(this->size(), Zero);
}


template<>
Foam::tmp<Foam::Field<Foam::label>>
Foam::transformFvPatchField<Foam::label>::gradientBoundaryCoeffs() const
{
    // TBD: Treat like zero-gradient
    return tmp<Field<label>>::New(this->size(), Zero);
}


// ************************************************************************* //
