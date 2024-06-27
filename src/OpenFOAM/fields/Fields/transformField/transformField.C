/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

#include "transformField.H"
#include "FieldM.H"
#include "diagTensor.H"

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

void Foam::transform
(
    vectorField& result,
    const quaternion& q,
    const vectorField& fld
)
{
    const tensor rot = q.R();
    if (result.cdata() == fld.cdata())
    {
        // std::for_each
        TSEQ_FORALL_F_OP_FUNC_S_F_inplace(result, =, transform, rot, fld)
    }
    else
    {
        // std::transform
        TSEQ_FORALL_F_OP_FUNC_S_F(result, =, transform, rot, fld)
    }
}


Foam::tmp<Foam::vectorField> Foam::transform
(
    const quaternion& q,
    const vectorField& fld
)
{
    auto tresult = tmp<vectorField>::New(fld.size());
    transform(tresult.ref(), q, fld);
    return tresult;
}


Foam::tmp<Foam::vectorField> Foam::transform
(
    const quaternion& q,
    const tmp<vectorField>& tfld
)
{
    tmp<vectorField> tresult = New(tfld);
    transform(tresult.ref(), q, tfld());
    tfld.clear();
    return tresult;
}


void Foam::transformPoints
(
    vectorField& result,
    const septernion& tr,
    const vectorField& fld
)
{
    const vector trans = tr.t();

    // Check if any translation
    if (mag(trans) > VSMALL)
    {
        if (result.cdata() == fld.cdata())
        {
            // std::for_each
            TSEQ_FORALL_F_OP_F_OP_S_inplace(result, =, fld, -, trans)
        }
        else
        {
            // std::transform
            TSEQ_FORALL_F_OP_F_OP_S(result, =, fld, -, trans)
        }
    }
    else
    {
        result = fld;
    }

    // Check if any rotation
    if (mag(tr.r().R() - I) > SMALL)
    {
        transform(result, tr.r(), result);
    }
}


Foam::tmp<Foam::vectorField> Foam::transformPoints
(
    const septernion& tr,
    const vectorField& fld
)
{
    auto tresult = tmp<vectorField>::New(fld.size());
    transformPoints(tresult.ref(), tr, fld);
    return tresult;
}


Foam::tmp<Foam::vectorField> Foam::transformPoints
(
    const septernion& tr,
    const tmp<vectorField>& tfld
)
{
    tmp<vectorField> tresult = New(tfld);
    transformPoints(tresult.ref(), tr, tfld());
    tfld.clear();
    return tresult;
}


// ************************************************************************* //
