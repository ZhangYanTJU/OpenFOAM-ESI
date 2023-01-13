/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "coupledFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makePatchFieldTypeNames(coupled);
    makePatchFieldTypeName(label, coupled);
}


// * * * * * * * * * * * * * * * Specialisations * * * * * * * * * * * * * * //

template<>
Foam::tmp<Foam::Field<Foam::label>>
Foam::coupledFvPatchField<Foam::label>::snGrad
(
    const scalarField& deltaCoeffs
) const
{
    // TBD: Treat like zero-gradient
    return tmp<Field<label>>::New(this->size(), Zero);
}


template<>
void Foam::coupledFvPatchField<Foam::label>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // TBD: Treat like zero-gradient
    fvPatchField<label>::operator=(this->patchInternalField());
    fvPatchField<label>::evaluate();
}


template<>
Foam::tmp<Foam::Field<Foam::label>>
Foam::coupledFvPatchField<Foam::label>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    // TBD: Treat like zero-gradient
    return tmp<Field<label>>::New(this->size(), label(1));
}


template<>
Foam::tmp<Foam::Field<Foam::label>>
Foam::coupledFvPatchField<Foam::label>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    // TBD: Treat like zero-gradient
    return tmp<Field<label>>::New(this->size(), Zero);
}


template<>
Foam::tmp<Foam::Field<Foam::label>>
Foam::coupledFvPatchField<Foam::label>::gradientInternalCoeffs
(
    const scalarField&
) const
{
    // TBD: Treat like zero-gradient
    return tmp<Field<label>>::New(this->size(), Zero);
}


template<>
Foam::tmp<Foam::Field<Foam::label>>
Foam::coupledFvPatchField<Foam::label>::gradientInternalCoeffs() const
{
    // TBD: Treat like zero-gradient
    return tmp<Field<label>>::New(this->size(), Zero);
}


// ************************************************************************* //
