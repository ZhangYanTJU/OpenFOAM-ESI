/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "cyclicFvPatchField.H"
#include "transformField.H"
#include "fvMatrixAssembly.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    cyclicPatch_(refCast<const cyclicFvPatch>(p))
{}


template<class Type>
Foam::cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    coupledFvPatchField<Type>(p, iF, dict, false), // Pass no valueRequired
    cyclicPatch_(refCast<const cyclicFvPatch>(p, dict))
{
    if (!isA<cyclicFvPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }

    if (valueRequired)
    {
        this->evaluate(Pstream::commsTypes::blocking);
    }
}


template<class Type>
Foam::cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const cyclicFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    cyclicPatch_(refCast<const cyclicFvPatch>(p))
{
    if (!isA<cyclicFvPatch>(this->patch()))
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
Foam::cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const cyclicFvPatchField<Type>& ptf
)
:
    cyclicLduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    cyclicPatch_(ptf.cyclicPatch_)
{}


template<class Type>
Foam::cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const cyclicFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(ptf, iF),
    cyclicPatch_(ptf.cyclicPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicFvPatchField<Type>::patchNeighbourField() const
{
    const Field<Type>& iField = this->primitiveField();
    const labelUList& nbrFaceCells =
        cyclicPatch().cyclicPatch().neighbPatch().faceCells();

    tmp<Field<Type>> tpnf(new Field<Type>(this->size()));
    Field<Type>& pnf = tpnf.ref();


    if (doTransform())
    {
        forAll(pnf, facei)
        {
            pnf[facei] = transform
            (
                forwardT()[0], iField[nbrFaceCells[facei]]
            );
        }
    }
    else
    {
        forAll(pnf, facei)
        {
            pnf[facei] = iField[nbrFaceCells[facei]];
        }
    }

    return tpnf;
}


template<class Type>
const Foam::cyclicFvPatchField<Type>&
Foam::cyclicFvPatchField<Type>::neighbourPatchField() const
{
    const GeometricField<Type, fvPatchField, volMesh>& fld =
    static_cast<const GeometricField<Type, fvPatchField, volMesh>&>
    (
        this->primitiveField()
    );

    return refCast<const cyclicFvPatchField<Type>>
    (
        fld.boundaryField()[this->cyclicPatch().neighbPatchID()]
    );
}



template<class Type>
void Foam::cyclicFvPatchField<Type>::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    const labelUList& nbrFaceCells =
        lduAddr.patchAddr
        (
            this->cyclicPatch().neighbPatchID()
        );

    solveScalarField pnf(psiInternal, nbrFaceCells);

    // Transform according to the transformation tensors
    transformCoupleField(pnf, cmpt);


    const labelUList& faceCells = lduAddr.patchAddr(patchId);

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, faceCells, coeffs, pnf);
}


template<class Type>
void Foam::cyclicFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes
) const
{
    const labelList& nbrFaceCells =
        lduAddr.patchAddr
        (
            this->cyclicPatch().neighbPatchID()
        );

    Field<Type> pnf(psiInternal, nbrFaceCells);

    // Transform according to the transformation tensors
    transformCoupleField(pnf);

    const labelUList& faceCells = lduAddr.patchAddr(patchId);

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, faceCells, coeffs, pnf);
}


template<class Type>
void Foam::cyclicFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
}


template<class Type>
void Foam::cyclicFvPatchField<Type>::manipulateMatrix
(
    fvMatrixAssembly& matrix,
    const labelList& faceMap,
    const label cellOffset
)
{
    const scalarField pAlphaSfDelta(gammaSfDelta());

    const labelUList& u = matrix.lduAddr().upperAddr();
    const labelUList& l = matrix.lduAddr().lowerAddr();

    DebugVar("manipulateMatrix")
    forAll (*this, faceI)
    {
        if (faceMap.size() > 0)
        {
            label globalFaceI = faceMap[faceI];
            if (globalFaceI != -1)
            {
                const scalar corr(pAlphaSfDelta[faceI]);
                if (cyclicPatch_.owner())
                {
                    matrix.lower()[globalFaceI] += corr;
                    matrix.upper()[globalFaceI] += corr;
                    matrix.diag()[u[globalFaceI]] -= corr;
                    matrix.diag()[l[globalFaceI]] -= corr;
                }
            }
            else
            {
                FatalErrorInFunction
                    << "Can not find faceId : " <<  globalFaceI
                    << exit(FatalError);
            }
        }
    }
}


template<class Type>
Foam::tmp<Foam::scalarField>
Foam::cyclicFvPatchField<Type>::gammaSfDelta() const
{
    const volScalarField& gamma =
        this->db().objectRegistry::template lookupObject
        <volScalarField>("(1|A(U))");

    const cyclicFvPatch& nbrPatch = cyclicPatch_.neighbPatch();

    scalarField gammab
    (
        gamma,
        cyclicPatch_.faceCells()
    );

    scalarField gammaNbr
    (
        gamma,
        nbrPatch.faceCells()
    );

    scalarField deltaf(1/cyclicPatch_.deltaCoeffs() + 1/nbrPatch.deltaCoeffs());

    scalarField gammaf
    (
        (
            gammab*(1/cyclicPatch_.deltaCoeffs())
          + gammaNbr*(1/nbrPatch.deltaCoeffs())
        )
        / deltaf
    );

    return
    (
         gammaf*this->patch().magSf()/deltaf
    );
}



// ************************************************************************* //
