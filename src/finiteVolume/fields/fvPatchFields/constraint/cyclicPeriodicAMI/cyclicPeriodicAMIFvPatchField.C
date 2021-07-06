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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::cyclicPeriodicAMIFvPatchField<Type>::cyclicPeriodicAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicPeriodicAMILduInterfaceField(),
    coupledFvPatchField<Type>(p, iF),
    cyclicPeriodicAMIPatch_(refCast<const cyclicPeriodicAMIFvPatch>(p))
{}


template<class Type>
Foam::cyclicPeriodicAMIFvPatchField<Type>::cyclicPeriodicAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    cyclicPeriodicAMILduInterfaceField(),
    coupledFvPatchField<Type>(p, iF, dict, dict.found("value")),
    cyclicPeriodicAMIPatch_(refCast<const cyclicPeriodicAMIFvPatch>(p, dict))
{
    if (!isA<cyclicPeriodicAMIFvPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }

    if (!dict.found("value"))
    {
        if (this->coupled())
        {
            this->evaluate(Pstream::commsTypes::blocking);
        }
        else
        {
            fvPatchField<Type>::operator=(this->patchInternalField());
        }
    }
}


template<class Type>
Foam::cyclicPeriodicAMIFvPatchField<Type>::cyclicPeriodicAMIFvPatchField
(
    const cyclicPeriodicAMIFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    cyclicPeriodicAMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    cyclicPeriodicAMIPatch_(refCast<const cyclicPeriodicAMIFvPatch>(p))
{
    if (!isA<cyclicPeriodicAMIFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::cyclicPeriodicAMIFvPatchField<Type>::cyclicPeriodicAMIFvPatchField
(
    const cyclicPeriodicAMIFvPatchField<Type>& ptf
)
:
    cyclicPeriodicAMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    cyclicPeriodicAMIPatch_(ptf.cyclicPeriodicAMIPatch_)
{}


template<class Type>
Foam::cyclicPeriodicAMIFvPatchField<Type>::cyclicPeriodicAMIFvPatchField
(
    const cyclicPeriodicAMIFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicPeriodicAMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf, iF),
    cyclicPeriodicAMIPatch_(ptf.cyclicPeriodicAMIPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//template<class Type>
//bool Foam::cyclicPeriodicAMIFvPatchField<Type>::coupled() const
//{
//    return cyclicPeriodicAMIPatch_.coupled();
//}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicPeriodicAMIFvPatchField<Type>::patchNeighbourField() const
{
    const Field<Type>& iField = this->primitiveField();
    const labelUList& nbrFaceCells =
        cyclicPeriodicAMIPatch_.neighbPatch().faceCells();

    Field<Type> pnf(iField, nbrFaceCells);

    if (debug&2)
    {
        Pout<< "cyclicPeriodicAMIFvPatchField<Type>::patchNeighbourField() :"
            << " field:" << this->internalField().name()
            << " patch:" << cyclicPeriodicAMIPatch_.name()
            << endl;
    }

    tmp<Field<Type>> tpnf;
    if (cyclicPeriodicAMIPatch_.applyLowWeightCorrection())
    {
        tpnf = cyclicPeriodicAMIPatch_.interpolate
        (
            pnf,
            this->patchInternalField()()
        );
    }
    else
    {
        tpnf = cyclicPeriodicAMIPatch_.interpolate(pnf);
    }

    if (this->doTransform())
    {
        tpnf.ref() = transform(this->forwardT(), tpnf());
    }
    return tpnf;
}


//template<class Type>
//const Foam::cyclicPeriodicAMIFvPatchField<Type>&
//Foam::cyclicPeriodicAMIFvPatchField<Type>::neighbourPatchField() const
//{
//    const GeometricField<Type, fvPatchField, volMesh>& fld =
//        static_cast<const GeometricField<Type, fvPatchField, volMesh>&>
//        (
//            this->primitiveField()
//        );
//
//    return refCast<const cyclicPeriodicAMIFvPatchField<Type>>
//    (
//        fld.boundaryField()[cyclicPeriodicAMIPatch_.neighbPatchID()]
//    );
//}


template<class Type>
void Foam::cyclicPeriodicAMIFvPatchField<Type>::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    // Problem: for vectorFields this routine gets called one component at
    // a time. Hence we cannot apply any transformation into cylindrical
    // coordinate system before we do the interpolation. Instead use the
    // component-coupled linear solver. In fvSolution:
    //    type            coupled;
    //    solver          PBiCCCG;    //PBiCGStab;

    const labelUList& nbrFaceCells =
        cyclicPeriodicAMIPatch_.neighbPatch().faceCells();

    solveScalarField pnf(psiInternal, nbrFaceCells);

    if (debug&2)
    {
        Pout<< "cyclicPeriodicAMIFvPatchField<Type>::updateInterfaceMatrix() :"
            << " field:" << this->internalField().name()
            << " patch:" << cyclicPeriodicAMIPatch_.name()
            << " cmpt:" << cmpt
            << endl;
    }

    // Transform according to the transformation tensors
    this->transformCoupleField(pnf, cmpt);

    if (cyclicPeriodicAMIPatch_.applyLowWeightCorrection())
    {
        solveScalarField pif(psiInternal, cyclicPeriodicAMIPatch_.faceCells());
        pnf = cyclicPeriodicAMIPatch_.interpolate(pnf, pif);
    }
    else
    {
        pnf = cyclicPeriodicAMIPatch_.interpolate(pnf);
    }

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, coeffs, pnf);
}


template<class Type>
void Foam::cyclicPeriodicAMIFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const bool add,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes
) const
{
    // For use in the component-coupled linear solver. In fvSolution:
    //    type            coupled;
    //    solver          PBiCCCG;    //PBiCGStab;

    const labelUList& nbrFaceCells =
        cyclicPeriodicAMIPatch_.neighbPatch().faceCells();

    Field<Type> pnf(psiInternal, nbrFaceCells);

    if (debug&2)
    {
        Pout<< "cyclicPeriodicAMIFvPatchField<Type>::updateInterfaceMatrix() :"
            << " field:" << this->internalField().name()
            << " patch:" << cyclicPeriodicAMIPatch_.name()
            << endl;
    }

    // Transform according to the transformation tensors
    this->transformCoupleField(pnf);

    if (cyclicPeriodicAMIPatch_.applyLowWeightCorrection())
    {
        Field<Type> pif(psiInternal, cyclicPeriodicAMIPatch_.faceCells());
        pnf = cyclicPeriodicAMIPatch_.interpolate(pnf, pif);
    }
    else
    {
        pnf = cyclicPeriodicAMIPatch_.interpolate(pnf);
    }

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, coeffs, pnf);
}


// ************************************************************************* //
