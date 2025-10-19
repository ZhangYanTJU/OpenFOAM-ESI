/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenCFD Ltd.
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
    Foam::HelicalForce

Group
    grpLagrangianIntermediateForceSubModels

Description
    This class implements a helical force model for Lagrangian particles.
    The force model applies a combination of centripetal, tangential, axial,
    and radial damping forces to the particles to induce helical motion around
    a specified axis.

Usage
    Minimal example by using \c constant/<cloud>/subModels/particleForces:
    \verbatim
    helical
    {
        // Mandatory entries
        axisPoint      <Function1<vector>>;  // (0 0 0);
        axisDir        <Function1<vector>>;  // (0 0 1);
        omega          <Function1<scalar>>;  // 10;   [rad/s]

        // Optional entries
        Uaxial         <Function1<scalar>>;  // 1;    [m/s]
        tauRotational  <Function1<scalar>>;  // 0.1;  [s]
        tauAxial       <Function1<scalar>>;  // 0.1;  [s]
        radialDamping  <Function1<scalar>>;  // 50;   [1/s]

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                  | Type | Reqd  | Deflt
      axisPoint | A point on the rotation axis | F1\<vector\> | yes | -
      axisDir   | Rotation axis (normalised)   | F1\<vector\> | yes | -
      omega | Rotational speed of the parcel   | F1\<scalar\> | yes | -
      Uaxial | Parcel axial speed in the rotation axis | F1\<scalar\> | no  | -
      tauRotational| Time-scale constant for tangential-speed feedback <!--
                --> | F1\<scalar\> | no | -
      tauAxial | Time-scale constant for axial-speed feedback  <!--
                --> | F1\<scalar\> | no | -
      radialDamping | Radial damping coefficient | F1\<scalar\> | no | -
    \endtable

Note
    - The \c axisDir vector must not be zero.
    - If \c tauRotational is specified, it must be positive.
    - If \c tauAxial is specified, it must be positive.
    - If \c tauAxial is specified, \c Uaxial must be specified.
    - If \c radialDamping is specified, it must be positive.

SourceFiles
    HelicalForce.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_HelicalForce_H
#define Foam_HelicalForce_H

#include "ParticleForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Type> class Function1;

/*---------------------------------------------------------------------------*\
                       Class HelicalForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class HelicalForce
:
    public ParticleForce<CloudType>
{
    // Private Data

        //- Point on the rotation axis, r0
        autoPtr<Function1<vector>> axisPointPtr_;

        //- Rotation axis (normalised), e
        autoPtr<Function1<vector>> axisDirPtr_;

        //- Rotational speed of the parcel about the rotation axis [rad/s]
        autoPtr<Function1<scalar>> omegaPtr_;

        //- Axial speed of the parcel in the rotation axis direction [m/s]
        autoPtr<Function1<scalar>> UaxialPtr_;

        //- Time-scale constant for the tangential feedback acceleration to
        //- keep the parcel's tangential velocity equal to the prescribed
        //- rotational speed around the axis, if need be [s]
        autoPtr<Function1<scalar>> tauRotationalPtr_;

        //- Time-scale constant for the axial feedback acceleration to keep
        //- the parcel's velocity along the axis equal to the prescribed axial
        //- speed, if need be [s]
        autoPtr<Function1<scalar>> tauAxialPtr_;

        //- Radial damping coefficient to suppress radial oscillations [1/s]
        //  e.g. 50 for strong damping
        autoPtr<Function1<scalar>> radialDampingPtr_;


public:

    //- Runtime type information
    TypeName("helical");


    // Constructors

        //- Construct from cloud, mesh and dictionary
        HelicalForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Copy construct
        HelicalForce(const HelicalForce<CloudType>& df);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new HelicalForce<CloudType>(*this)
            );
        }

        //- No copy assignment
        void operator=(const HelicalForce<CloudType>&) = delete;


    //- Destructor
    virtual ~HelicalForce();


    // Member Functions

    // Evaluation

        //- Calculate the coupled force
        virtual forceSuSp calcCoupled
        (
            const typename CloudType::parcelType& p,
            const typename CloudType::parcelType::trackingData& td,
            const scalar dt,
            const scalar mass,
            const scalar Re,
            const scalar muc
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HelicalForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
