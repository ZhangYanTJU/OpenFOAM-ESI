/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "externalFileSource.H"
#include "faMatrices.H"
#include "faCFD.H"
#include "zeroGradientFaPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fa
{
    defineTypeNameAndDebug(externalFileSource, 0);

    addToRunTimeSelectionTable
    (
        option,
        externalFileSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fa::externalFileSource::externalFileSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvPatch& p
)
:
    faceSetOption(sourceName, modelType, dict, p),
    fieldName_(dict.get<word>("fieldName")),
    tableName_(dict.get<word>("tableName")),
    value_
    (
        new PatchFunction1Types::MappedFile<scalar>
        (
            p.patch(),
            "uniformValue",
            dict,
            tableName_,          // field table name
            true                 // face values
        )
    ),
    curTimeIndex_(-1)
{
    fieldNames_.setSize(1, fieldName_);

    applied_.setSize(fieldNames_.size(), false);
  
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fa::externalFileSource::~externalFileSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fa::externalFileSource::addSup
(
    const areaScalarField& solidMass,
    faMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (isActive())
    {
        DebugInfo<< name() << ": applying source to " << eqn.psi().name()<<endl;
        
        const scalar t = mesh().time().timeOutputValue();
        
        if (curTimeIndex_ != mesh().time().timeIndex())
        {
            IOobject io
            (
                "Q",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            );

            areaScalarField p
            (
                io,
                regionMesh(),
                dimensionedScalar("p", dimPressure, Zero),
                zeroGradientFaPatchScalarField::typeName
            );
            
            p.field() = value_->value(t);
            eqn += p/solidMass;
            curTimeIndex_ = mesh().time().timeIndex();
        }
    }
}


bool Foam::fa::externalFileSource::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        return true;
    }

    return false;
}

// ************************************************************************* //

