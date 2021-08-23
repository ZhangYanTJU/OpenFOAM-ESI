/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "volFields.H"
#include "cellCellStencil.H"
#include "cellCellStencilObject.H"
#include "oversetFvMeshBase.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::oversetFvPatchField<Type>::oversetFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    zeroGradientFvPatchField<Type>(p, iF),
    oversetPatch_(refCast<const oversetFvPatch>(p)),
    setHoleCellValue_(false),
    interpolateHoleCellValue_(false),
    holeCellValue_(pTraits<Type>::min)
{}


template<class Type>
Foam::oversetFvPatchField<Type>::oversetFvPatchField
(
    const oversetFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchField<Type>(ptf, p, iF, mapper),
    oversetPatch_(refCast<const oversetFvPatch>(p)),
    setHoleCellValue_(ptf.setHoleCellValue_),
    interpolateHoleCellValue_(ptf.interpolateHoleCellValue_),
    holeCellValue_(ptf.holeCellValue_)
{}


template<class Type>
Foam::oversetFvPatchField<Type>::oversetFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchField<Type>(p, iF, dict),
    oversetPatch_(refCast<const oversetFvPatch>(p, dict)),
    setHoleCellValue_(dict.lookupOrDefault<bool>("setHoleCellValue", false)),
    interpolateHoleCellValue_
    (
        dict.lookupOrDefault<bool>("interpolateHoleCellValue", false)
    ),
    holeCellValue_
    (
        setHoleCellValue_
      ? dict.get<Type>("holeCellValue")
      : pTraits<Type>::min
    )
{}


template<class Type>
Foam::oversetFvPatchField<Type>::oversetFvPatchField
(
    const oversetFvPatchField<Type>& ptf
)
:
    zeroGradientFvPatchField<Type>(ptf),
    oversetPatch_(ptf.oversetPatch_),
    setHoleCellValue_(ptf.setHoleCellValue_),
    interpolateHoleCellValue_(ptf.interpolateHoleCellValue_),
    holeCellValue_(ptf.holeCellValue_)
{}


template<class Type>
Foam::oversetFvPatchField<Type>::oversetFvPatchField
(
    const oversetFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    zeroGradientFvPatchField<Type>(ptf, iF),
    oversetPatch_(ptf.oversetPatch_),
    setHoleCellValue_(ptf.setHoleCellValue_),
    interpolateHoleCellValue_(ptf.interpolateHoleCellValue_),
    holeCellValue_(ptf.holeCellValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::oversetFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (this->oversetPatch_.master())
    {
        // Trigger interpolation
        const fvMesh& mesh = this->internalField().mesh();
        const dictionary& fvSchemes = mesh.schemesDict();
        const word& fldName = this->internalField().name();

        if (&mesh.lduAddr() != &mesh.fvMesh::lduAddr())
        {
            // Running extended addressing. Flux correction already done
            // in the linear solver (through the initUpdateInterfaceMatrix
            // method below)
            if (debug)
            {
                Info<< "Skipping overset interpolation for solved-for field "
                    << fldName << endl;
            }
        }
        else if (!fvSchemes.found("oversetInterpolation"))
        {
            IOWarningInFunction(fvSchemes)
                << "Missing required dictionary entry"
                << " 'oversetInterpolation'"
                << ". Skipping overset interpolation for field "
                << fldName << endl;
        }
        else if (fvSchemes.found("oversetInterpolationRequired"))
        {
            // Backwards compatibility mode: only interpolate what is
            // explicitly mentioned

            if (fvSchemes.found("oversetInterpolationSuppressed"))
            {
                FatalIOErrorInFunction(fvSchemes)
                    << "Cannot have both dictionary entry"
                    << " 'oversetInterpolationSuppresed' and "
                    << " 'oversetInterpolationRequired' for field "
                    << fldName << exit(FatalIOError);
            }
            const dictionary& intDict = fvSchemes.subDict
            (
                "oversetInterpolationRequired"
            );
            if (intDict.found(fldName))
            {
                if (debug)
                {
                    Info<< "Interpolating field " << fldName << endl;
                }

                // Interpolate without bc update (since would trigger infinite
                // recursion back to oversetFvPatchField<Type>::evaluate)
                // The only problem is bcs that use the new cell values
                // (e.g. zeroGradient, processor). These need to appear -after-
                // the 'overset' bc.
                mesh.interpolate
                (
                    const_cast<Field<Type>&>
                    (
                        this->primitiveField()
                    )
                );
            }
            else if (debug)
            {
                Info<< "Skipping overset interpolation for field "
                    << fldName << endl;
            }
        }
        else
        {
            const dictionary* dictPtr
            (
                fvSchemes.findDict("oversetInterpolationSuppressed")
            );

            const wordHashSet& suppress =
                Stencil::New(mesh).nonInterpolatedFields();

            bool skipInterpolate = suppress.found(fldName);

            if (dictPtr)
            {
                skipInterpolate =
                    skipInterpolate
                 || dictPtr->found(fldName);
            }

            if (skipInterpolate)
            {
                if (debug)
                {
                    Info<< "Skipping suppressed overset interpolation"
                        << " for field " << fldName << endl;
                }
            }
            else
            {
                if (debug)
                {
                    Info<< "Interpolating non-suppressed field " << fldName
                        << endl;
                }

                // Interpolate without bc update (since would trigger infinite
                // recursion back to oversetFvPatchField<Type>::evaluate)
                // The only problem is bcs that use the new cell values
                // (e.g. zeroGradient, processor). These need to appear -after-
                // the 'overset' bc.

                const cellCellStencil& overlap = Stencil::New(mesh);

                Field<Type>& fld = const_cast<Field<Type>&>(this->primitiveField());

                // tbd: different weights for different variables ...
                cellCellStencil::interpolate
                (
                    fld,
                    mesh,
                    overlap,
                    overlap.cellInterpolationWeights()
                );

                if (this->setHoleCellValue_)
                {
                    const labelUList& types = overlap.cellTypes();
                    label nConstrained = 0;
                    forAll(types, celli)
                    {
                        const label cType = types[celli];
                        if
                        (
                            cType == cellCellStencil::HOLE
                        || cType == cellCellStencil::SPECIAL
                        )
                        {
                            fld[celli] = this->holeCellValue_;
                            nConstrained++;
                        }
                    }

                    if (debug)
                    {
                        Pout<< FUNCTION_NAME << " field:" << fldName
                            << " patch:" << this->oversetPatch_.name()
                            << " set:" << nConstrained << " cells to:"
                            << this->holeCellValue_ << endl;
                    }
                }
                /*
                mesh.interpolate
                (
                    const_cast<Field<Type>&>
                    (
                        this->primitiveField()
                    )
                );
                */
            }
        }
    }

    zeroGradientFvPatchField<Type>::initEvaluate(commsType);
}


template<class Type>
void Foam::oversetFvPatchField<Type>::manipulateMatrix
(
    fvMatrix<Type>& matrix
)
{
    //const word& fldName = this->internalField().name();

    if (this->manipulatedMatrix())
    {
        return;
    }

    const oversetFvPatch& ovp = this->oversetPatch_;

    if (ovp.master())
    {
        // Trigger interpolation
        const fvMesh& mesh = this->internalField().mesh();
        const word& fldName = this->internalField().name();

        // Try to find out if the solve routine comes from the unadapted mesh
        // TBD. This should be cleaner.
        if (&mesh.lduAddr() == &mesh.fvMesh::lduAddr())
        {
            //Pout<< FUNCTION_NAME << "SKIPPING MANIP field:" << fldName
            //    << " patch:" << ovp.name() << endl;
            return;
        }

        if (debug)
        {
            Pout<< FUNCTION_NAME << " field:" << fldName
                << " patch:" << ovp.name() << endl;
        }


        // Calculate stabilised diagonal as normalisation for interpolation
        const scalarField norm
        (
            dynamic_cast<const oversetFvMeshBase&>(mesh).normalisation(matrix)
        );

        const cellCellStencil& overlap = Stencil::New(mesh);
        const labelUList& types = overlap.cellTypes();
        const labelListList& stencil = overlap.cellStencil();

        dynamic_cast<const oversetFvMeshBase&>(mesh).addInterpolation
        (
            matrix,
            norm,
            this->setHoleCellValue_,
            this->holeCellValue_
        );

        if (debug)
        {
            pointField allCs(mesh.cellCentres());
            const mapDistribute& map = overlap.cellInterpolationMap();
            map.distribute(allCs, false, UPstream::msgType()+1);

            // Make sure we don't interpolate from a hole

            scalarField marker(this->primitiveField().size(), 0);
            forAll(types, celli)
            {
                if (types[celli] == cellCellStencil::HOLE)
                {
                    marker[celli] = 1.0;
                }
            }
            cellCellStencil::interpolate
            (
                marker,
                mesh,
                overlap,
                overlap.cellInterpolationWeights()
            );

            forAll(marker, celli)
            {
                if
                (
                    types[celli] == cellCellStencil::INTERPOLATED
                 && marker[celli] > SMALL
                )
                {
                    //FatalErrorInFunction
                    WarningInFunction
                        << " field:" << fldName
                        << " patch:" << ovp.name()
                        << " found:" << celli
                        << " at:" << mesh.cellCentres()[celli]
                        << " donorSlots:" << stencil[celli]
                        << " at:"
                        << UIndirectList<point>(allCs, stencil[celli])
                        << " amount-of-hole:" << marker[celli]
                        //<< exit(FatalError);
                        << endl;
                }
            }

            // Make sure we don't have matrix coefficients for interpolated
            // or hole cells

            const lduAddressing& addr = mesh.lduAddr();
            const labelUList& upperAddr = addr.upperAddr();
            const labelUList& lowerAddr = addr.lowerAddr();
            const scalarField& lower = matrix.lower();
            const scalarField& upper = matrix.upper();

            forAll(lowerAddr, facei)
            {
                const label l = lowerAddr[facei];
                const bool lHole = (types[l] == cellCellStencil::HOLE);
                const label u = upperAddr[facei];
                const bool uHole = (types[u] == cellCellStencil::HOLE);

                if
                (
                    (lHole && upper[facei] != 0.0)
                 || (uHole && lower[facei] != 0.0)
                )
                {
                    FatalErrorInFunction
                        << "Hole-neighbouring face:" << facei
                        << " lower:" << l
                        << " type:" << types[l]
                        << " coeff:" << lower[facei]
                        << " upper:" << upperAddr[facei]
                        << " type:" << types[u]
                        << " coeff:" << upper[facei]
                        << exit(FatalError);
                }


                // Empty donor list: treat like hole but still allow to
                // influence neighbouring domains
                const bool lEmpty =
                (
                    types[l] == cellCellStencil::INTERPOLATED
                 && stencil[l].empty()
                );
                const bool uEmpty =
                (
                    types[u] == cellCellStencil::INTERPOLATED
                 && stencil[u].empty()
                );

                if
                (
                    (lEmpty && upper[facei] != 0.0)
                 || (uEmpty && lower[facei] != 0.0)
                )
                {
                    FatalErrorInFunction
                        << "Still connected face:" << facei << " lower:" << l
                        << " type:" << types[l]
                        << " coeff:" << lower[facei]
                        << " upper:" << u
                        << " type:" << types[u]
                        << " coeff:" << upper[facei]
                        << exit(FatalError);
                }
            }

            forAll(matrix.internalCoeffs(), patchi)
            {
                const labelUList& fc = addr.patchAddr(patchi);
                //const Field<Type>& intCoeffs =
                //    matrix.internalCoeffs()[patchi];
                const Field<Type>& bouCoeffs = matrix.boundaryCoeffs()[patchi];
                forAll(fc, i)
                {
                    label celli = fc[i];

                    const bool lHole = (types[celli] == cellCellStencil::HOLE);
                    if (lHole && bouCoeffs[i] != pTraits<Type>::zero)
                    {
                        FatalErrorInFunction
                            << "Patch:" << patchi
                            << " patchFace:" << i
                            << " lower:" << celli
                            << " type:" << types[celli]
                            << " bouCoeff:" << bouCoeffs[i]
                            << exit(FatalError);
                    }

                    // Check whether I am influenced by neighbouring domains
                    const bool lEmpty =
                    (
                        types[celli] == cellCellStencil::INTERPOLATED
                     && stencil[celli].empty()
                    );

                    if (lEmpty && bouCoeffs[i] != pTraits<Type>::zero)
                    {
                        FatalErrorInFunction
                            << "Patch:" << patchi
                            << " patchFace:" << i
                            << " lower:" << celli
                            << " type:" << types[celli]
                            << " bouCoeff:" << bouCoeffs[i]
                            << exit(FatalError);
                    }
                }
            }


            // Make sure that diagonal is non-zero. Note: should add
            // boundaryCoeff ...
            const FieldField<Field, Type>& internalCoeffs =
                matrix.internalCoeffs();
            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
            {
                // Replacement for m.addBoundaryDiag(norm, cmpt);
                scalarField diag(matrix.diag());
                forAll(internalCoeffs, patchi)
                {
                    const labelUList& fc = addr.patchAddr(patchi);
                    const Field<Type>& intCoeffs = internalCoeffs[patchi];
                    const scalarField cmptCoeffs(intCoeffs.component(cmpt));
                    forAll(fc, i)
                    {
                        diag[fc[i]] += cmptCoeffs[i];
                    }
                }

                forAll(diag, celli)
                {
                    if (mag(diag[celli]) < SMALL)
                    {
                        FatalErrorInFunction
                            << "Patch:" << ovp.name()
                            << " cell:" << celli
                            << " at:" << mesh.cellCentres()[celli]
                            << " diag:" << diag[celli]
                            << exit(FatalError);
                    }
                }
            }
        }
    }

    zeroGradientFvPatchField<Type>::manipulateMatrix(matrix);
}


template<class Type>
void Foam::oversetFvPatchField<Type>::write(Ostream& os) const
{
    zeroGradientFvPatchField<Type>::write(os);
    if (this->setHoleCellValue_)
    {
        os.writeEntry("setHoleCellValue", setHoleCellValue_);
        os.writeEntry("holeCellValue", holeCellValue_);
        os.writeEntryIfDifferent
        (
            "interpolateHoleCellValue",
            false,
            interpolateHoleCellValue_
        );
    }
    // Make sure to write the value for ease of postprocessing.
    this->writeEntry("value", os);
}


// ************************************************************************* //
