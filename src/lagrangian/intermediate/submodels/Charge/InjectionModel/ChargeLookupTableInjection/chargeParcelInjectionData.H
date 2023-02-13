/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

Class
    Foam::chargeParcelInjectionData

Description
    Container class to provide injection data for electric charge parcels

SourceFiles
    chargeParcelInjectionData.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_chargeParcelInjectionData_H
#define Foam_chargeParcelInjectionData_H

#include "reactingMultiphaseParcelInjectionData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class chargeParcelInjectionData;

// Forward declaration of friend functions

Ostream& operator<<
(
    Ostream&,
    const chargeParcelInjectionData&
);

Istream& operator>>
(
    Istream&,
    chargeParcelInjectionData&
);

/*---------------------------------------------------------------------------*\
            Class chargeParcelInjectionData Declaration
\*---------------------------------------------------------------------------*/

class chargeParcelInjectionData
:
    public reactingMultiphaseParcelInjectionData
{
protected:

    // Parcel properties

        //- Electric charge [C]
        scalar eq_;


public:

    //- Runtime type information
    TypeName("chargeParcelInjectionData");

    // Constructors

        //- Null constructor
        chargeParcelInjectionData();

        //- Construct from dictionary
        chargeParcelInjectionData(const dictionary& dict);

        //- Construct from Istream
        chargeParcelInjectionData(Istream& is);


    //- Destructor
    virtual ~chargeParcelInjectionData();


    // Access

        //- Return const access to the electric charge
        inline scalar eq() const;


    // Edit

        //- Return access to the electric charge
        inline scalar& eq();


    // I-O

        //- Ostream operator
        friend Ostream& operator<<
        (
            Ostream& os,
            const chargeParcelInjectionData& data
        );

        //- Istream operator
        friend Istream& operator>>
        (
            Istream& is,
            chargeParcelInjectionData& data
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "chargeParcelInjectionDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
