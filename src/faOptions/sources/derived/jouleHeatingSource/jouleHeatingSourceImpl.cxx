/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2025 OpenCFD Ltd.
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

#include "emptyFaPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::fa::jouleHeatingSource::initialiseSigma
(
    const dictionary& dict,
    autoPtr<Function1<Type>>& sigmaFunctionPtr
)
{
    typedef GeometricField<Type, faPatchField, areaMesh> FieldType;

    const word sigmaName
    (
        IOobject::scopedName(typeName, "sigma") + suffixHint()
    );

    auto& obr = regionMesh().thisDb();

    IOobject io
    (
        sigmaName,
        obr.time().timeName(),
        obr,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        IOobject::REGISTER
    );

    autoPtr<FieldType> tsigma;

    // Is sigma defined using a Function1 type?
    sigmaFunctionPtr = Function1<Type>::NewIfPresent("sigma", dict, &mesh_);

    if (sigmaFunctionPtr)
    {
        tsigma.reset
        (
            new FieldType
            (
                io,
                regionMesh(),
                Foam::zero{},  // value
                sqr(dimCurrent)/dimPower/dimLength
            )
        );

        Info<< "    Conductivity 'sigma' read from dictionary as f(T)"
            << nl << endl;
    }
    else
    {
        // Sigma to be defined by user input
        io.readOpt(IOobject::MUST_READ);

        tsigma.reset(new FieldType(io, regionMesh()));

        Info<< "    Conductivity 'sigma' read from file" << nl << endl;
    }

    regIOobject::store(tsigma);
}


template<class Type>
const Foam::GeometricField<Type, Foam::faPatchField, Foam::areaMesh>&
Foam::fa::jouleHeatingSource::updateSigma
(
    const autoPtr<Function1<Type>>& sigmaFunctionPtr
) const
{
    typedef GeometricField<Type, faPatchField, areaMesh> FieldType;

    const word sigmaName
    (
        IOobject::scopedName(typeName, "sigma") + suffixHint()
    );

    const auto& obr = regionMesh().thisDb();

    auto& sigma = obr.lookupObjectRef<FieldType>(sigmaName);

    if (!sigmaFunctionPtr)
    {
        // Electrical conductivity field, sigma, was specified by the user
        return sigma;
    }

    const auto& sigmaFunction = sigmaFunctionPtr();

    const auto& T = obr.lookupObject<areaScalarField>(TName_);

    // Internal field
    forAll(sigma, i)
    {
        sigma[i] = sigmaFunction.value(T[i]);
    }


    // Boundary field
    auto& bf = sigma.boundaryFieldRef();
    forAll(bf, patchi)
    {
        faPatchField<Type>& pf = bf[patchi];
        if (!isA<emptyFaPatch>(pf))
        {
            const scalarField& Tbf = T.boundaryField()[patchi];
            forAll(pf, facei)
            {
                pf[facei] = sigmaFunction.value(Tbf[facei]);
            }
        }
    }

    // Update processor patches
    sigma.correctBoundaryConditions();

    return sigma;
}


// ************************************************************************* //
