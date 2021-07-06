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

\*---------------------------------------------------------------------------*/

#include "cyclicPeriodicAMIGAMGInterfaceField.H"
#include "addToRunTimeSelectionTable.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicPeriodicAMIGAMGInterfaceField, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        cyclicPeriodicAMIGAMGInterfaceField,
        lduInterface
    );
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        cyclicPeriodicAMIGAMGInterfaceField,
        lduInterfaceField
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cyclicPeriodicAMIGAMGInterfaceField::cyclicPeriodicAMIGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const lduInterfaceField& fineInterface
)
:
    GAMGInterfaceField(GAMGCp, fineInterface),
    cyclicPeriodicAMIInterface_
    (
        refCast<const cyclicPeriodicAMIGAMGInterface>(GAMGCp)
    ),
    doTransform_(false),
    rank_(0)
{
    const cyclicPeriodicAMILduInterfaceField& p =
        refCast<const cyclicPeriodicAMILduInterfaceField>(fineInterface);

    doTransform_ = p.doTransform();
    rank_ = p.rank();
}


Foam::cyclicPeriodicAMIGAMGInterfaceField::cyclicPeriodicAMIGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const bool doTransform,
    const int rank
)
:
    GAMGInterfaceField(GAMGCp, doTransform, rank),
    cyclicPeriodicAMIInterface_
    (
        refCast<const cyclicPeriodicAMIGAMGInterface>(GAMGCp)
    ),
    doTransform_(doTransform),
    rank_(rank)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicPeriodicAMIGAMGInterfaceField::
~cyclicPeriodicAMIGAMGInterfaceField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cyclicPeriodicAMIGAMGInterfaceField::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    const auto& neighbPatch = cyclicPeriodicAMIInterface_.neighbPatch();

    // Get neighbouring field
    solveScalarField pnf
    (
        neighbPatch.interfaceInternalField
        (
            psiInternal
        )
    );

    // Transform according to the transformation tensors
    transformCoupleField(pnf, cmpt);

    if (debug&2)
    {
        Pout<< "cyclicPeriodicAMIGAMGInterfaceField::updateInterfaceMatrix :"
            //<< " field:" << this->internalField().name()
            //<< " interface:" << cyclicPeriodicAMIInterface_.name()
            << " cmpt:" << cmpt
            << endl;
    }

    if (cyclicPeriodicAMIInterface_.owner())
    {
        pnf = cyclicPeriodicAMIInterface_.AMI().interpolateToSource(pnf);
    }
    else
    {
        pnf = neighbPatch.AMI().interpolateToTarget(pnf);
    }

    this->addToInternalField(result, !add, coeffs, pnf);
}


// ************************************************************************* //
