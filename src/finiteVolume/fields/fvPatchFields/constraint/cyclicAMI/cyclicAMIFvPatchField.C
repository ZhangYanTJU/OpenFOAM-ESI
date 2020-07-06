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
        const scalarField pAlphaSfDelta(gammaSfDelta(matrix, iMatrix));

        const labelUList& u = matrix.lduAddr().upperAddr();
        const labelUList& l = matrix.lduAddr().lowerAddr();

        label subFaceI = 0;

        if (faceMap.size() > 0)
        {
            forAll (faceMap, j)
            {
                label globalFaceI = faceMap[j];

                if (globalFaceI != -1)
                {
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
            }
        }
    }
}

template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::manipulateInterBoundCoeffs
(
    fvMatrixAssembly& matrix,
    const scalarField& weights,
    const labelUList& fc,
    const vectorField& delta,
    const label iMatrix,
    const label neigProc,
    scalarField& boundaryCoeffs,
    scalarField& internalCoeffs
)
{
    const volScalarField& gamma =
        this->db().objectRegistry::template lookupObject
        <
            volScalarField
        >(diffusivityName_);

    const scalarListList& srcWeight =
        cyclicAMIPatch_.cyclicAMIPatch().AMI().srcWeights();

    const labelListList& sourceFaces =
        cyclicAMIPatch_.cyclicAMIPatch().AMI().srcAddress();

    const cyclicAMIFvPatch& neigPatch = cyclicAMIPatch_.neighbPatch();

    const tmpNrc<mapDistribute> tgtMapPtr =
        cyclicAMIPatch_.cyclicAMIPatch().AMI().tgtMap();

    List<scalarField> srcFaceToTgtGamma;
    collectStencilData
    (
        tgtMapPtr,
        sourceFaces,
        scalarField(gamma, neigPatch.faceCells()),
        srcFaceToTgtGamma
    );
    const scalarField& magSf = this->patch().magSf();

    scalarField gammaf(fc.size(), Zero);
    scalarField Sf(fc.size(), Zero);

    const labelList& rmtFaces =
        matrix.mesh().myRmtTgtFaces()[iMatrix][this->patch().index()][neigProc];

    const labelList& faceCellsIds = cyclicAMIPatch_.faceCells();

    forAll (rmtFaces, subFaceI)
    {
        label faceI = rmtFaces[subFaceI];

        forAll(sourceFaces[faceI], j)
        {
            label cellI = faceCellsIds[faceI];

            const scalarList& w = srcWeight[faceI];

            {
                //Remote
                gammaf[subFaceI] =
                    gamma[cellI]*(1 - weights[subFaceI])
                  + weights[subFaceI]*srcFaceToTgtGamma[faceI][j];

                Sf[subFaceI] = w[j]*magSf[faceI];

            }
        }
    }

    boundaryCoeffs = -gammaf*Sf/(mag(delta) + ROOTSMALL);
    internalCoeffs = -gammaf*Sf/(mag(delta) + ROOTSMALL);
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

    const labelList& nbrCellIds =
        matrix.mesh().cellBoundMap()[iMatrix][index];
    const labelList& cellIds =
        matrix.mesh().cellBoundMap()[iMatrix][nbrPathID];

    const labelList& faceMap =
        matrix.mesh().faceBoundMap()[iMatrix][index];

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
            const label localFaceId =
                matrix.mesh().magSfFaceBoundMap()[iMatrix][index][subFaceI];

            if (localFaceId != -1)
            {
                Sf[subFaceI] = w[i]*this->patch().magSf()[localFaceId];

                vector unitArea = patchSf[localFaceId]/magSf[localFaceId];

                //nonOrthDeltaCoeffs from uncorrected SnGrad
                deltaf[subFaceI] =
                    max
                    (
                        unitArea & deltafV[subFaceI],
                        0.05*mag(deltafV[subFaceI])
                    );
            }
            subFaceI++;
        }
    }

    const scalarField gammaSfDelta(gammaf*Sf/deltaf);

    // Set internalCoeffs and boundaryCoeffs in the assembly matrix
    // on clyclicAMI patches to be used in the individual matrix by
    // matrix.flux()

    if (matrix.mesh().fluxRequired(this->internalField().name()))
    {
        const label globalPatchID =
            matrix.mesh().patchLocalToGlobalMap()[iMatrix][index];

        matrix.internalCoeffs().set(globalPatchID, -gammaSfDelta);
        matrix.boundaryCoeffs().set(globalPatchID, -gammaSfDelta);

        const label nbrGlobalPatchID =
            matrix.mesh().patchLocalToGlobalMap()[iMatrix][nbrPathID];

        matrix.internalCoeffs().set(nbrGlobalPatchID, -gammaSfDelta);
        matrix.boundaryCoeffs().set(nbrGlobalPatchID, -gammaSfDelta);
    }

    return tmp<scalarField>(new scalarField(gammaSfDelta));
}


template<class Type>
template<class Type2>
void Foam::cyclicAMIFvPatchField<Type>::collectStencilData
(
    const tmpNrc<mapDistribute>& mapPtr,
    const labelListList& stencil,
    const Type2& data,
    List<Type2>& expandedData
)
{
    expandedData.setSize(stencil.size());
    if (mapPtr.valid())
    {
        Type2 work(data);
        mapPtr().distribute(work);

        forAll(stencil, facei)
        {
            const labelList& slots = stencil[facei];
            expandedData[facei].append
            (
                UIndirectList<typename Type2::value_type>(work, slots)
            );
        }
    }
    else
    {
        forAll(stencil, facei)
        {
            const labelList& slots = stencil[facei];
            expandedData[facei].append
            (
                UIndirectList<typename Type2::value_type>(data, slots)
            );
        }
    }
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeEntryIfDifferent<word>("diffusivity", "none", diffusivityName_);
    this->writeEntry("value", os);
}


// ************************************************************************* //

