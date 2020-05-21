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

#include "fvMatrixAssembly.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(p, iF),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p)),
    diffusivityName_("none")
{}


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(p, iF, dict, dict.found("value")),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p, dict)),
    diffusivityName_(dict.lookupOrDefault<word>("diffusivity", "none"))
{
    if (!isA<cyclicAMIFvPatch>(p))
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
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const cyclicAMIFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p)),
    diffusivityName_(ptf.diffusivityName_)
{
    if (!isA<cyclicAMIFvPatch>(this->patch()))
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
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const cyclicAMIFvPatchField<Type>& ptf
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    cyclicAMIPatch_(ptf.cyclicAMIPatch_),
    diffusivityName_(ptf.diffusivityName_)
{}


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const cyclicAMIFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf, iF),
    cyclicAMIPatch_(ptf.cyclicAMIPatch_),
    diffusivityName_(ptf.diffusivityName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::cyclicAMIFvPatchField<Type>::coupled() const
{
    return cyclicAMIPatch_.coupled();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicAMIFvPatchField<Type>::patchNeighbourField() const
{
    const Field<Type>& iField = this->primitiveField();
    const labelUList& nbrFaceCells =
        cyclicAMIPatch_.cyclicAMIPatch().neighbPatch().faceCells();

    Field<Type> pnf(iField, nbrFaceCells);

    tmp<Field<Type>> tpnf;
    if (cyclicAMIPatch_.applyLowWeightCorrection())
    {
        tpnf = cyclicAMIPatch_.interpolate(pnf, this->patchInternalField()());
    }
    else
    {
        tpnf = cyclicAMIPatch_.interpolate(pnf);
    }

    if (doTransform())
    {
        tpnf.ref() = transform(forwardT(), tpnf());
    }

    return tpnf;
}


template<class Type>
const Foam::cyclicAMIFvPatchField<Type>&
Foam::cyclicAMIFvPatchField<Type>::neighbourPatchField() const
{
    const GeometricField<Type, fvPatchField, volMesh>& fld =
        static_cast<const GeometricField<Type, fvPatchField, volMesh>&>
        (
            this->primitiveField()
        );

    return refCast<const cyclicAMIFvPatchField<Type>>
    (
        fld.boundaryField()[cyclicAMIPatch_.neighbPatchID()]
    );
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    const labelUList& nbrFaceCells =
        lduAddr.patchAddr
        (
            cyclicAMIPatch_.cyclicAMIPatch().neighbPatchID()
        );

    solveScalarField pnf(psiInternal, nbrFaceCells);

    // Transform according to the transformation tensors
    transformCoupleField(pnf, cmpt);

    if (cyclicAMIPatch_.applyLowWeightCorrection())
    {
        solveScalarField pif(psiInternal, cyclicAMIPatch_.faceCells());
        pnf = cyclicAMIPatch_.interpolate(pnf, pif);
    }
    else
    {
        pnf = cyclicAMIPatch_.interpolate(pnf);
    }

    const labelUList& faceCells = lduAddr.patchAddr(patchId);

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, faceCells, coeffs, pnf);
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::updateInterfaceMatrix
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
    const labelUList& nbrFaceCells =
        lduAddr.patchAddr
        (
            cyclicAMIPatch_.cyclicAMIPatch().neighbPatchID()
        );

    Field<Type> pnf(psiInternal, nbrFaceCells);

    // Transform according to the transformation tensors
    transformCoupleField(pnf);

    if (cyclicAMIPatch_.applyLowWeightCorrection())
    {
        Field<Type> pif(psiInternal, cyclicAMIPatch_.faceCells());
        pnf = cyclicAMIPatch_.interpolate(pnf, pif);
    }
    else
    {
        pnf = cyclicAMIPatch_.interpolate(pnf);
    }

    const labelUList& faceCells = lduAddr.patchAddr(patchId);

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, faceCells, coeffs, pnf);
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::manipulateMatrix
(
    fvMatrixAssembly& matrix,
    const labelList& faceMap,
    const label cellOffset,
    const label iMatrix
)
{
    if (cyclicAMIPatch_.owner())
    {
        /*
        const volScalarField& divPhiByA =
            this->db().objectRegistry::template lookupObject
            <volScalarField>("divPhiByA");

        const volScalarField& divPhiByAold = divPhiByA.oldTime();

        const label index(this->patch().index());
        label nbrPathID = cyclicAMIPatch_.cyclicAMIPatch().neighbPatchID();

        const labelList& cellIds =
            matrix.mesh().cellBoundMap()[iMatrix][index];

        const labelList& nbrCellIds =
            matrix.mesh().cellBoundMap()[iMatrix][nbrPathID];

        tmp<scalarField> V(matrix.mesh().meshes()[iMatrix].V());

        scalar relax = 1;
        forAll (cellIds, i)
        {
            label cellI = cellIds[i];
            matrix.source()[cellI] =
                relax*matrix.source()[cellI] + (1-relax)*divPhiByAold[cellI]*V()[cellI];
        }
        forAll (nbrCellIds, i)
        {
            label nbrCellI = nbrCellIds[i];
            matrix.source()[nbrCellI] =
                relax*matrix.source()[nbrCellI]
                + (1-relax)*divPhiByAold[nbrCellI]*V()[nbrCellI];
        }
        */

        //const scalarListList& srcWeight =
        //    cyclicAMIPatch_.cyclicAMIPatch().AMI().srcWeights();

        const scalarField pAlphaSfDelta(gammaSfDelta(matrix, iMatrix));

        const labelUList& u = matrix.lduAddr().upperAddr();
        const labelUList& l = matrix.lduAddr().lowerAddr();

        label subFaceI(0);
        //forAll (*this, faceI)
        {
            if (faceMap.size() > 0)
            {
                //const scalarList& w = srcWeight[faceI];

                //for(label i=0; i<w.size(); i++)
                forAll (faceMap, j)
                {
                    label globalFaceI = faceMap[j];
//DebugVar(globalFaceI)
                    if (globalFaceI != -1)
                    {
//DebugVar(subFaceI)
                        const scalar corr(pAlphaSfDelta[subFaceI]);

                        matrix.upper()[globalFaceI] += corr;
                        matrix.diag()[u[globalFaceI]] -= corr;
                        matrix.diag()[l[globalFaceI]] -= corr;

                        if (matrix.asymmetric())
                        {
                            matrix.lower()[globalFaceI] += corr;
                        }
                        subFaceI++;
                    }
//                     else
//                     {
//                         FatalErrorInFunction
//                             << "Can not find faceId : " <<  globalFaceI
//                             << exit(FatalError);
//                     }

//                     subFaceI++;
                }
            }
        }
    }
}

template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::manipulateInterBoundCoeffs
(
    fvMatrixAssembly& matrix,
    const labelList& faceMap,
    const label cellOffset,
    const label iMatrix,
    scalarField& boundaryCoeffs,
    scalarField& internalCoeffs
)
{
    DebugVar("manipulateInterBoundCoeffs")
}


template<class Type>
Foam::tmp<Foam::scalarField>
Foam::cyclicAMIFvPatchField<Type>::gammaSfDelta
(
    fvMatrixAssembly& matrix,
    const label iMatrix
) const
{
    const label index(this->patch().index());

    const volScalarField& gamma =
        this->db().objectRegistry::template lookupObject
        <
            volScalarField
        >(diffusivityName_);

    label nbrPathID = cyclicAMIPatch_.cyclicAMIPatch().neighbPatchID();

    const labelList& nbrCellIds = matrix.mesh().cellBoundMap()[iMatrix][index];
    const labelList& cellIds = matrix.mesh().cellBoundMap()[iMatrix][nbrPathID];

    const labelList& faceMap = matrix.mesh().faceBoundMap()[iMatrix][index];
    label nSubFaces = 0;
    forAll (faceMap, i)
    {
        if (faceMap[i] != -1)
        {
            nSubFaces++;
        }
    }

    scalarField gammaf(nSubFaces, Zero);
    scalarField deltaf(nSubFaces, Zero);
    vectorField deltafV(nSubFaces, Zero);
    scalarField Sf(nSubFaces, Zero);

    tmp<vectorField> C(matrix.mesh().meshes()[iMatrix].C());

    label subFaceI = 0;
    forAll (cellIds, i)
    {
        label cellI = cellIds[i];
        label nbrCellI = nbrCellIds[i];

        if (cellI != -1 && nbrCellI != -1)
        {
            gammaf[subFaceI] = 0.5*(gamma[cellI] + gamma[nbrCellI]);

            // Outgoing from the patch
            deltafV[subFaceI] = C()[nbrCellI] - C()[cellI];
        }
        subFaceI++;
    }

    const scalarListList& srcWeight =
        cyclicAMIPatch_.cyclicAMIPatch().AMI().srcWeights();

    const vectorField& patchSf = this->patch().Sf();
    const scalarField& magSf = this->patch().magSf();

    subFaceI = 0;
    forAll(*this, faceI)
    {
        const scalarList& w = srcWeight[faceI];

        for(label i=0; i<w.size(); i++)
        {
            const label faceIdMap =
                matrix.mesh().magSfFaceBoundMap()[iMatrix][index][subFaceI];
            if (faceIdMap != -1)
            {
                Sf[subFaceI] = w[i]*this->patch().magSf()[faceIdMap];

                vector unitArea = patchSf[faceIdMap]/magSf[faceIdMap];

                //nonOrthDeltaCoeffs from corrected SnGrad
                deltaf[subFaceI] =
                    max(unitArea & deltafV[subFaceI], 0.05*mag(deltafV[subFaceI]));

                //NonOrthCorrectionVectors
                //corrVecs[subFaceI] =
                //    unitArea - deltafV[subFaceI]/deltaf[subFaceI];
            }

            subFaceI++;
        }
    }

    scalarField gammaSfDelta(gammaf*Sf/deltaf);

//DebugVar(gammaSfDelta);

    // Set internalCoeffs and boundaryCoeffs in the assembly matrix
    // on clyclicAMI patches to be used in the individual matrix by
    // matrix.flux()

    if (matrix.mesh().fluxRequired(this->internalField().name()))
    {
        const label globalPatchID =
            matrix.mesh().patchLocalToGlobalMap()[iMatrix][index];
//DebugVar(index)
//DebugVar(globalPatchID)
        matrix.internalCoeffs().set(globalPatchID, -gammaSfDelta);
        matrix.boundaryCoeffs().set(globalPatchID, -gammaSfDelta);

        const label nbrGlobalPatchID =
            matrix.mesh().patchLocalToGlobalMap()[iMatrix][nbrPathID];
//DebugVar(nbrPathID)
//DebugVar(nbrGlobalPatchID)
        matrix.internalCoeffs().set(nbrGlobalPatchID, -gammaSfDelta);
        matrix.boundaryCoeffs().set(nbrGlobalPatchID, -gammaSfDelta);
    }

    return tmp<scalarField>(new scalarField(gammaSfDelta));
}

template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeEntry("diffusivity", diffusivityName_);
    this->writeEntry("value", os);
}


// ************************************************************************* //

