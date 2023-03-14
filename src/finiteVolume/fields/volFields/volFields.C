/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTemplateTypeNameAndDebug(volScalarField::Internal, 0);
defineTemplateTypeNameAndDebug(volVectorField::Internal, 0);
defineTemplateTypeNameAndDebug(volSphericalTensorField::Internal, 0);
defineTemplateTypeNameAndDebug(volSymmTensorField::Internal, 0);
defineTemplateTypeNameAndDebug(volTensorField::Internal, 0);

defineTemplateTypeNameAndDebug(volScalarField, 0);
defineTemplateTypeNameAndDebug(volVectorField, 0);
defineTemplateTypeNameAndDebug(volSphericalTensorField, 0);
defineTemplateTypeNameAndDebug(volSymmTensorField, 0);
defineTemplateTypeNameAndDebug(volTensorField, 0);

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Specializations

namespace Foam
{

template<>
tmp<GeometricField<scalar, fvPatchField, volMesh>>
GeometricField<scalar, fvPatchField, volMesh>::component
(
    const direction
) const
{
    return *this;
}


template<>
void GeometricField<scalar, fvPatchField, volMesh>::replace
(
    const direction,
    const GeometricField<scalar, fvPatchField, volMesh>& gsf
)
{
    *this == gsf;
}


template<>
bool GeometricBoundaryField<scalar, fvPatchField, volMesh>::check
(
    const scalar tol,
    const bool doExit
) const
{
    if (!this->size())
    {
        return true;
    }

    if (GeometricField<scalar, fvPatchField, volMesh>::debug)
    {
        const auto& fvp0 = this->operator[](0);
        PoutInFunction
            << "Checking boundary consistency for field "
            << fvp0.internalField().name()
            << endl;
    }

    // Store old value
    List<Field<scalar>> oldBfld(this->size());
    boolList oldUpdated(this->size());
    boolList oldManipulated(this->size());

    forAll(*this, patchi)
    {
        if (this->set(patchi) && this->operator[](patchi).size())
        {
            const auto& fvp = this->operator[](patchi);
            oldUpdated[patchi] = fvp.updated();
            oldManipulated[patchi] = fvp.manipulatedMatrix();
            oldBfld[patchi] = fvp;
        }
    }

    // Re-evaluate
    const_cast<GeometricBoundaryField<scalar, fvPatchField, volMesh>&>
    (
        *this
    ).evaluate();

    // Check
    bool ok = true;
    forAll(*this, patchi)
    {
        if (this->set(patchi) && this->operator[](patchi).size())
        {
            const auto& pfld = this->operator[](patchi);
            const auto& oldPfld = oldBfld[patchi];

            //if (pfld != oldBfld[patchi])
            forAll(pfld, facei)
            {
                if (mag(pfld[facei]-oldPfld[facei]) > tol)
                {
                    ok = false;

                    if (doExit)
                    {
                        FatalErrorInFunction << "Field "
                            << pfld.internalField().name()
                            << " is not evaluated?"
                            << " On patch " << pfld.patch().name()
                            << " : average of field = "
                            << average(oldBfld[patchi])
                            << ". Average of evaluated field = "
                            << average(pfld)
                            << exit(FatalError);
                    }
                }
            }
        }
    }

    // Restore bfld, updated
    forAll(*this, patchi)
    {
        if (this->set(patchi) && this->operator[](patchi).size())
        {
            auto& fvp = const_cast<fvPatchField<scalar>&>
            (
                this->operator[](patchi)
            );
            fvp.setUpdated(oldUpdated[patchi]);
            fvp.setManipulated(oldManipulated[patchi]);
            Field<scalar>& vals = fvp;
            vals = std::move(oldBfld[patchi]);
        }
    }

    if (GeometricField<scalar, fvPatchField, volMesh>::debug)
    {
        const auto& fvp0 = this->operator[](0);
        PoutInFunction
            << "Result of checking for field "
            << fvp0.internalField().name()
            << " : " << ok
            << endl;
    }

    return ok;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Global Data * * * * * * * * * * * * * * * //

// Note hard-coded values are more reliable than other alternatives

const Foam::wordList Foam::fieldTypes::internal
({
    "volScalarField::Internal",
    "volVectorField::Internal",
    "volSphericalTensorField::Internal",
    "volSymmTensorField::Internal",
    "volTensorField::Internal"
});


const Foam::wordList Foam::fieldTypes::volume
({
    "volScalarField",
    "volVectorField",
    "volSphericalTensorField",
    "volSymmTensorField",
    "volTensorField"
});


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
