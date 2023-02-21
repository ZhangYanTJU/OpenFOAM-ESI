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
    Foam::SpaceChargeDensity

Group
    grpLagrangianIntermediateFunctionObjects

Description

    References:
    \verbatim
        Governing equations (tag:YSSD):
            Ye, Q., Steigleder, T., Scheibe, A., & Domnick, J. (2002).
            Numerical simulation of the electrostatic powder coating process
            with a corona spray gun. Journal of Electrostatics, 54(2), 189-205.
            DOI:10.1016/S0304-3886(01)00181-4
    \endverbatim

SourceFiles
    SpaceChargeDensity.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_SpaceChargeDensity_H
#define Foam_SpaceChargeDensity_H

#include "CloudFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class SpaceChargeDensity Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SpaceChargeDensity
:
    public CloudFunctionObject<CloudType>
{
    // Private Data

        // Typedefs

            //- Convenience typedef for parcel type
            typedef typename CloudType::parcelType parcelType;

        //- Particle mass field
        autoPtr<volScalarField> mPtr_;

        //- Particle electric-charge field
        autoPtr<volScalarField> qPtr_;

        //- Particle specific concentration field
        autoPtr<volScalarField> cPtr_;

        //- Space charge density field
        autoPtr<volScalarField> erhoPtr_;


protected:

    // Protected Member Functions

        //- Write post-processing info
        virtual void write();


public:

    //- Runtime type information
    TypeName("spaceChargeDensity");


    // Constructors

        //- Construct from dictionary
        SpaceChargeDensity
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Copy construct
        SpaceChargeDensity(const SpaceChargeDensity<CloudType>& vf);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new SpaceChargeDensity<CloudType>(*this)
            );
        }

        //- No copy assignment
        void operator=(const SpaceChargeDensity<CloudType>&) = delete;


    //- Destructor
    virtual ~SpaceChargeDensity() = default;


    // Member Functions

        // Evaluation

        //- Pre-evolve hook
        virtual void preEvolve
        (
            const typename parcelType::trackingData& td
        );

        //- Post-evolve hook
        virtual void postEvolve
        (
            const typename parcelType::trackingData& td
        );

        //- Post-move hook
        virtual void postMove
        (
            parcelType& p,
            const scalar dt,
            const point& position0,
            bool& keepParticle
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SpaceChargeDensity.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
