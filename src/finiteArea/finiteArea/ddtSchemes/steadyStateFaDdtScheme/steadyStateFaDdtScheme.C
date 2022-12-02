/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

#include "steadyStateFaDdtScheme.H"
#include "facDiv.H"
#include "faMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
steadyStateFaDdtScheme<Type>::facDdt
(
    const dimensioned<Type> dt
)
{
    return tmp<GeometricField<Type, faPatchField, areaMesh>>::New
    (
        this->fieldIOobject("ddt("+dt.name()+')'),
        mesh(),
        dimensioned<Type>(dt.dimensions()/dimTime, Zero)
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
steadyStateFaDdtScheme<Type>::facDdt0
(
    const dimensioned<Type> dt
)
{
    return tmp<GeometricField<Type, faPatchField, areaMesh>>::New
    (
        this->fieldIOobject("ddt("+dt.name()+')'),
        mesh(),
        dimensioned<Type>(dt.dimensions()/dimTime, Zero)
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
steadyStateFaDdtScheme<Type>::facDdt
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return tmp<GeometricField<Type, faPatchField, areaMesh>>::New
    (
        this->fieldIOobject("ddt("+vf.name()+')'),
        mesh(),
        dimensioned<Type>(vf.dimensions()/dimTime, Zero)
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
steadyStateFaDdtScheme<Type>::facDdt0
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return tmp<GeometricField<Type, faPatchField, areaMesh>>::New
    (
        this->fieldIOobject("ddt0("+vf.name()+')'),
        mesh(),
        dimensioned<Type>(vf.dimensions()/dimTime, Zero)
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
steadyStateFaDdtScheme<Type>::facDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return tmp<GeometricField<Type, faPatchField, areaMesh>>::New
    (
        this->fieldIOobject("ddt("+rho.name()+','+vf.name()+')'),
        mesh(),
        dimensioned<Type>(rho.dimensions()*vf.dimensions()/dimTime, Zero)
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
steadyStateFaDdtScheme<Type>::facDdt0
(
    const dimensionedScalar& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return tmp<GeometricField<Type, faPatchField, areaMesh>>::New
    (
        this->fieldIOobject("ddt0("+rho.name()+','+vf.name()+')'),
        mesh(),
        dimensioned<Type>(rho.dimensions()*vf.dimensions()/dimTime, Zero)
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
steadyStateFaDdtScheme<Type>::facDdt
(
    const areaScalarField& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return tmp<GeometricField<Type, faPatchField, areaMesh>>::New
    (
        this->fieldIOobject("ddt("+rho.name()+','+vf.name()+')'),
        mesh(),
        dimensioned<Type>(rho.dimensions()*vf.dimensions()/dimTime, Zero)
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
steadyStateFaDdtScheme<Type>::facDdt0
(
    const areaScalarField& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return tmp<GeometricField<Type, faPatchField, areaMesh>>::New
    (
        this->fieldIOobject("ddt0("+rho.name()+','+vf.name()+')'),
        mesh(),
        dimensioned<Type>(rho.dimensions()*vf.dimensions()/dimTime, Zero)
    );
}


template<class Type>
tmp<faMatrix<Type>>
steadyStateFaDdtScheme<Type>::famDdt
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return tmp<faMatrix<Type>>::New
    (
        vf,
        vf.dimensions()*dimArea/dimTime
    );
}


template<class Type>
tmp<faMatrix<Type>>
steadyStateFaDdtScheme<Type>::famDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return tmp<faMatrix<Type>>::New
    (
        vf,
        rho.dimensions()*vf.dimensions()*dimArea/dimTime
    );
}


template<class Type>
tmp<faMatrix<Type>>
steadyStateFaDdtScheme<Type>::famDdt
(
    const areaScalarField& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return tmp<faMatrix<Type>>::New
    (
        vf,
        rho.dimensions()*vf.dimensions()*dimArea/dimTime
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
