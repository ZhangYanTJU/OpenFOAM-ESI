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

Class
    Foam::fa::jouleHeatingSource

Group
    grpFvOptionsSources

Description
    Evolves an electrical potential equation

    \f[
        \grad \left( \sigma \grad V \right)
    \f]

    where \f$ V \f$ is electrical potential and \f$\sigma\f$ is the
    electrical current

    To provide a Joule heating contribution according to:

    Differential form of Joule heating - power per unit volume:

    \f[
        \frac{d(P)}{d(V)} = J \cdot E
    \f]

    where \f$ J \f$ is the current density and \f$ E \f$ the electric field.
    If no magnetic field is present:

    \f[
        J = \sigma E
    \f]

    The electric field given by

    \f[
        E = \grad V
    \f]

    Therefore:

    \f[
        \frac{d(P)}{d(V)} = J \cdot E
                          = (sigma E) \cdot E
                          = (sigma \grad V) \cdot \grad V
    \f]


Usage
    Isotropic (scalar) electrical conductivity
    \verbatim
    jouleHeatingSourceCoeffs
    {
        heating
        {
            type                                jouleHeatingSource;
            anisotropicElectricalConductivity    no;
            T                                   Ts_ceilingShell;

            active                              true;
            region                              ceilingShell;
            selectionMode                       all;
        }
    }
    \endverbatim


    Where:
    \table
        Property     | Description               | Required  | Default value
        T            | Name of temperature field | yes
        sigma        | Electrical conductivity as a function of temperature |no|
        anisotropicElectricalConductivity | Anisotropic flag | yes |
        region       | 2D Region to be applied   | yes
    \endtable

    The electrical conductivity can be specified using either:
    - If the \c sigma entry is present the electrical conductivity is specified
      as a function of temperature using a Function1 type
    - If not present the sigma field will be read from file
    - If the anisotropicElectricalConductivity flag is set to 'true', sigma
      should be specified as a tensor quantity

SourceFiles
    jouleHeatingSource.C

SeeAlso
    Foam::Function1

\*---------------------------------------------------------------------------*/

#ifndef fa_jouleHeatingSource_H
#define fa_jouleHeatingSource_H

#include "faOption.H"
#include "Function1.H"
#include "areaFields.H"
#include "faceSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fa
{

/*---------------------------------------------------------------------------*\
                        Class jouleHeatingSource Declaration
\*---------------------------------------------------------------------------*/

class jouleHeatingSource
:
    public faceSetOption
{
    // Private data

        //- Name of temperature field - default = "T" (optional)
        word TName_;

        //- Electrical potential field / [V]
        areaScalarField V_;

        //- Flag to indicate that the electrical conductivity is anisotropic
        bool anisotropicElectricalConductivity_;

        //- Electrical conductivity as a scalar function of temperature
        autoPtr<Function1<scalar>> scalarSigmaVsTPtr_;

        //- Electrical conductivity as a vector function of temperature
        autoPtr<Function1<tensor>> tensorSigmaVsTPtr_;

        //- Current time index (used for updating)
        label curTimeIndex_;


    // Private Member Functions

        //- No copy construct
        jouleHeatingSource(const jouleHeatingSource&) = delete;

        //- No copy assignment
        void operator=(const jouleHeatingSource&) = delete;

        //- Initialise the electrical conductivity field
        template<class Type>
        void initialiseSigma
        (
            const dictionary& dict,
            autoPtr<Function1<Type>>& sigmaVsTPtr
        );

        //- Update the electrical conductivity field
        template<class Type>
        const GeometricField<Type, faPatchField, areaMesh>&
        updateSigma(const autoPtr<Function1<Type>>& sigmaVsTPtr) const;


public:

    //- Runtime type information
    TypeName("jouleHeatingSource");


    // Constructors

        //- Construct from explicit source name and mesh
        jouleHeatingSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh,
            const fvPatch& patch
        );


    //- Destructor
    virtual ~jouleHeatingSource();


    // Member Functions

        // Evaluate

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const areaScalarField& h,
                const areaScalarField& rho,
                faMatrix<scalar>& eqn,
                const label fieldi
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "jouleHeatingSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
