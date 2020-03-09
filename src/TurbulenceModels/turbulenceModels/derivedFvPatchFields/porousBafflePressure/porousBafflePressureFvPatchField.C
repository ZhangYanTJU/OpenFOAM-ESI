/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "porousBafflePressureFvPatchField.H"
#include "surfaceFields.H"
#include "turbulenceModel.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrixAssembly.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousBafflePressureFvPatchField::porousBafflePressureFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpFvPatchField<scalar>(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    D_(),
    I_(),
    length_(0),
    uniformJump_(false)
{}


Foam::porousBafflePressureFvPatchField::porousBafflePressureFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedJumpFvPatchField<scalar>(p, iF),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    rhoName_(dict.getOrDefault<word>("rho", "rho")),
    D_(Function1<scalar>::New("D", dict)),
    I_(Function1<scalar>::New("I", dict)),
    length_(dict.get<scalar>("length")),
    uniformJump_(dict.getOrDefault("uniformJump", false))
{
    fvPatchField<scalar>::operator=
    (
        Field<scalar>("value", dict, p.size())
    );
}


Foam::porousBafflePressureFvPatchField::porousBafflePressureFvPatchField
(
    const porousBafflePressureFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedJumpFvPatchField<scalar>(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    D_(ptf.D_.clone()),
    I_(ptf.I_.clone()),
    length_(ptf.length_),
    uniformJump_(ptf.uniformJump_)
{}


Foam::porousBafflePressureFvPatchField::porousBafflePressureFvPatchField
(
    const porousBafflePressureFvPatchField& ptf
)
:
    cyclicLduInterfaceField(),
    fixedJumpFvPatchField<scalar>(ptf),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    D_(ptf.D_.clone()),
    I_(ptf.I_.clone()),
    length_(ptf.length_),
    uniformJump_(ptf.uniformJump_)
{}


Foam::porousBafflePressureFvPatchField::porousBafflePressureFvPatchField
(
    const porousBafflePressureFvPatchField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpFvPatchField<scalar>(ptf, iF),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    D_(ptf.D_.clone()),
    I_(ptf.I_.clone()),
    length_(ptf.length_),
    uniformJump_(ptf.uniformJump_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porousBafflePressureFvPatchField::manipulateMatrix
(
    fvMatrixAssembly& matrix,
    const labelList& faceMap,
    const label cellOffset
)
{
    const cyclicPolyPatch& cpp =
        refCast<const cyclicPolyPatch>(patch().patch());

    //const scalarField vf(gammaDeltaVf());
    const scalarField pAlphaSfDelta(gammaSfDelta());

//DebugVar(pAlphaSfDelta)

    const labelUList& u = matrix.lduAddr().upperAddr();
    const labelUList& l = matrix.lduAddr().lowerAddr();

    forAll (*this, faceI)
    {
        if (faceMap.size() > 0)
        {
            label globalFaceI = faceMap[faceI];
            if (globalFaceI != -1)
            {
                //const scalar corr(vf[faceI]*pAlphaSfDelta[faceI]);
                const scalar corr(pAlphaSfDelta[faceI]);
                if (cpp.owner())
                {
/*                {
DebugVar(globalFaceI)
DebugVar(u[globalFaceI])
DebugVar(l[globalFaceI])
DebugVar(matrix.lower()[globalFaceI])
DebugVar(matrix.upper()[globalFaceI])
DebugVar(matrix.diag()[u[globalFaceI]])
DebugVar(matrix.diag()[l[globalFaceI]])
*/
                    matrix.lower()[globalFaceI] += corr;
                    matrix.upper()[globalFaceI] += corr;
                    matrix.diag()[u[globalFaceI]] -= corr;
                    matrix.diag()[l[globalFaceI]] -= corr;
/*
DebugVar(matrix.lower()[globalFaceI])
DebugVar(matrix.upper()[globalFaceI])
DebugVar(matrix.diag()[u[globalFaceI]])
DebugVar(matrix.diag()[l[globalFaceI]])
*/

                }
                /*
                else
                {
                    matrix.lower()[globalFaceI] -= corr;
                    matrix.diag()[u[globalFaceI]] += corr;
                }
                */
                /*
                if (cpp.owner())
                {
                    const labelUList& l = matrix.lduAddr().lowerAddr();
                    matrix.diag()[l[globalFaceI]] += corr;
                    matrix.lower()[globalFaceI] -= corr;
                }
                else
                {
                    matrix.upper()[globalFaceI] -= corr;
                    const labelUList& u = matrix.lduAddr().upperAddr();
                    matrix.diag()[u[globalFaceI]] += corr;
                }
                */
            }
            else
            {
                FatalErrorInFunction
                    << "Can not find faceId : " <<  globalFaceI
                    << exit(FatalError);
            }
        }
    }
    /*
    const scalarField sourceCorrection
    (
        pAlphaSfDelta
        *(
            deltaH()*vf
          + deltaQflux()/beta()
        )
    );

    const labelUList& fc = patch().faceCells();

    forAll(fc, i)
    {
        label localCelli = fc[i];
        label globalCelli = cellOffset + localCelli;
        matrix.source()[globalCelli] += sourceCorrection[i];
    }
    */
}

Foam::tmp<Foam::scalarField>
Foam::porousBafflePressureFvPatchField::gammaDeltaVf() const
{

    const label nbrPatchi = this->cyclicPatch().neighbPatchID();

    const label patchId = patch().index();

    //const label patchi = patch().index();
    //const polyMesh& nbrMesh = mpp.sampleMesh();
    const volScalarField& gamma =
            db().lookupObject<volScalarField>("(1|A(U))");

    const cyclicFvPatch& nbrPatch = this->cyclicPatch().neighbPatch();

    scalarField gammaDeltaNbr
    (
        gamma.boundaryField()[nbrPatchi]*nbrPatch.deltaCoeffs()
    );


    scalarField gammaDelta
    (
        gamma.boundaryField()[patchId]*patch().deltaCoeffs()
    );

    return (gammaDeltaNbr/(gammaDeltaNbr + gammaDelta));
}


Foam::tmp<Foam::scalarField>
Foam::porousBafflePressureFvPatchField::gammaSfDelta() const
{
    const volScalarField& gamma =
            db().lookupObject<volScalarField>("(1|A(U))");

    const cyclicFvPatch& nbrPatch = this->cyclicPatch().neighbPatch();

    scalarField gammab
    (
        gamma,
        this->cyclicPatch().faceCells()
    );
//DebugVar(this->cyclicPatch().faceCells())

    scalarField gammaNbr
    (
        gamma,
        nbrPatch.faceCells()
    );

    scalarField deltaf(1/patch().deltaCoeffs() + 1/nbrPatch.deltaCoeffs());

    scalarField gammaf
    (
        (gammab*(1/patch().deltaCoeffs()) + gammaNbr*(1/nbrPatch.deltaCoeffs()))
      /  deltaf
    );

//DebugVar(gammaf)

    return
    (
         gammaf*patch().magSf()/deltaf
    );
}



void Foam::porousBafflePressureFvPatchField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const surfaceScalarField& phi =
            db().lookupObject<surfaceScalarField>(phiName_);

    const fvsPatchField<scalar>& phip =
        patch().patchField<surfaceScalarField, scalar>(phi);

    scalarField Un(phip/patch().magSf());

    if (phi.dimensions() == dimMass/dimTime)
    {
        Un /= patch().lookupPatchField<volScalarField, scalar>(rhoName_);
    }

    if (uniformJump_)
    {
        Un = gAverage(Un);
    }
    scalarField magUn(mag(Un));

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const scalar t = db().time().timeOutputValue();
    const scalar D = D_->value(t);
    const scalar I = I_->value(t);

    setJump
    (
        -sign(Un)
        *(
            D*turbModel.nu(patch().index())
          + I*0.5*magUn
         )*magUn*length_
    );

    if (internalField().dimensions() == dimPressure)
    {
        setJump
        (
            jump()*patch().lookupPatchField<volScalarField, scalar>(rhoName_)
        );
    }

    if (debug)
    {
        scalar avePressureJump = gAverage(jump());
        scalar aveVelocity = gAverage(Un);

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << " Average pressure drop :" << avePressureJump
            << " Average velocity :" << aveVelocity
            << endl;
    }

    fixedJumpFvPatchField<scalar>::updateCoeffs();
}


void Foam::porousBafflePressureFvPatchField::write(Ostream& os) const
{
    fixedJumpFvPatchField<scalar>::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    D_->writeData(os);
    I_->writeData(os);
    os.writeEntry("length", length_);
    os.writeEntry("uniformJump", uniformJump_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        porousBafflePressureFvPatchField
    );
}

// ************************************************************************* //
