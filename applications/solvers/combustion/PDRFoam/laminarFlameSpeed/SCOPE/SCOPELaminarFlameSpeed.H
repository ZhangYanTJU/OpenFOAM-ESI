/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::laminarFlameSpeedModels::SCOPE

Description
    Laminar flame speed obtained from the SCOPE correlation.

    Seven parameters are specified in terms of polynomial functions of
    stoichiometry. Two polynomials are fitted, covering different parts of the
    flammable range. If the mixture is outside the fitted range, linear
    interpolation is used between the extreme of the polynomio and the upper or
    lower flammable limit with the Markstein number constant.

    Variations of pressure and temperature from the reference values are taken
    into account through \f$ pexp \f$ and \f$ texp \f$

    The laminar burning velocity fitting polynomial is:

    \f$ Su = a_{0}(1+a_{1}x+K+..a_{i}x^{i}..+a_{6}x^{6}) (p/p_{ref})^{pexp}
    (T/T_{ref})^{texp} \f$

    where:

        \f$ a_{i} \f$ are the polinomial coefficients.

        \f$ pexp \f$ and \f$ texp \f$ are the pressure and temperature factors
        respectively.

        \f$ x \f$ is the equivalence ratio.

        \f$ T_{ref} \f$ and \f$ p_{ref} \f$ are the temperature and pressure
        references for the laminar burning velocity.


SourceFiles
    SCOPELaminarFlameSpeed.C

\*---------------------------------------------------------------------------*/

#ifndef SCOPE_H
#define SCOPE_H

#include "laminarFlameSpeed.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarFlameSpeedModels
{

/*---------------------------------------------------------------------------*\
                           Class SCOPE Declaration
\*---------------------------------------------------------------------------*/

class SCOPE
:
    public laminarFlameSpeed
{
    // Private Data

        class polynomial
        :
            public FixedList<scalar, 7>
        {
        public:

            //- Lower limit
            scalar ll;

            //- Upper polynomial limit
            scalar ul;

            //- Value at lower limit
            scalar llv;

            //- Value at upper limit
            scalar ulv;

            //- Changeover point from lower to upper polynomial
            scalar lu;

            //- Construct from dictionary
            polynomial(const dictionary& polyDict);
        };


        dictionary coeffsDict_;

        //- Lower flammability limit
        scalar LFL_;

        //- Upper flammability limit
        scalar UFL_;

        //- Lower Su polynomial
        polynomial SuPolyL_;

        //- Upper Su polynomial
        polynomial SuPolyU_;

        //- Temperature correction exponent0
        scalar Texp0_;

        //- Temperature correction exponent1
        scalar Texp1_;

        //- Temperature correction exponent2
        scalar Texp2_;

        //- Temperature correction exponent3
        scalar Texp3_;

        //- Temperature correction exponent4
        scalar Texp4_;

        //- Temperature correction exponent5
        scalar Texp5_;

        //- Temperature correction exponent6
        scalar Texp6_;

        //- Pressure correction exponent0
        scalar pexp0_;

        //- Pressure correction exponent1
        scalar pexp1_;

        //- pressure correction exponent2
        scalar pexp2_;

        //- pressure correction exponent3
        scalar pexp3_;

        //- Pressure correction exponent4
        scalar pexp4_;

        //- Pressure correction exponent5
        scalar pexp5_;

        //- pressure correction exponent6
        scalar pexp6_;

        //- Additional coefficient
        scalar CIn_;

        //- Lower Ma polynomial
        polynomial MaPolyL_;

        //- Upper Ma polynomial
        polynomial MaPolyU_;


    // Private Member Functions

        //- Polynomial evaluated from the given equivalence ratio
        //- and polynomial coefficients
        static inline scalar polyPhi(scalar phi, const polynomial& a);

        //- Laminar flame speed evaluated from the given equivalence ratio
        //- at the reference temperature and pressure
        inline scalar SuRef(scalar phi) const;

        //- Markstein evaluated from the given equivalence ratio
        inline scalar Ma(scalar phi) const;

        //- Laminar flame speed evaluated from the given equivalence ratio
        //- corrected for temperature and pressure dependence
        inline scalar Su0pTphi(scalar p, scalar Tu, scalar phi) const;

        //- Laminar flame speed evaluated from the given uniform
        //- equivalence ratio corrected for temperature and pressure dependence
        tmp<volScalarField> Su0pTphi
        (
            const volScalarField& p,
            const volScalarField& Tu,
            scalar phi
        ) const;

        //- Laminar flame speed evaluated from the given equivalence ratio
        //- distribution corrected for temperature and pressure dependence
        tmp<volScalarField> Su0pTphi
        (
            const volScalarField& p,
            const volScalarField& Tu,
            const volScalarField& phi
        ) const;

        //- Return the Markstein number
        //- evaluated from the given equivalence ratio
        tmp<volScalarField> Ma(const volScalarField& phi) const;

        //- No copy construct
        SCOPE(const SCOPE&) = delete;

        //- No copy assign
        void operator=(const SCOPE&) = delete;


public:

    //- Runtime type information
    TypeName("SCOPE");

    // Constructors

        //- Construct from dictionary and psiuReactionThermo
        SCOPE
        (
            const dictionary& dict,
            const psiuReactionThermo& ct
        );


    //- Destructor
    ~SCOPE();


    // Member functions

        //- Return the Markstein number
        tmp<volScalarField> Ma() const;

        //- Return CIn
        scalar CIn() const;

        //- Return the laminar flame speed [m/s]
        tmp<volScalarField> operator()() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End laminarFlameSpeedModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
