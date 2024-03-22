/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::evapotranspirationHeatTransferModels::tree

Description
    Applies sources on temperature (\c T - incompressible) or energy
    (\c h/e - compressible) equation to incorporate evapotranspiration
    heat-transfer effects from the specified tree canopy.

    The heat-transfer is calculated based on the following relations.
    The heat transfer is given by the first equation and the
    evapotranspiration is calculated linearly proportional to absorbed
    solar radiation (direct and diffusive) with the empirical relation
    shown in the second equation below.

    \f[
        Q = \lambda E_t \frac{1000}{3600} \int_{cellZone} LAD \, dV
    \f]

    where
    \vartable
      Q         | Heat transfer from tree crown                 [W]
      \lambda   | Specific latent heat of vaporisation          [MJ/kg]
      E_t       | Area-averaged evapotranspiration              [g/m^2/hr]
      LAD       | Leaf area density                             [m^2/m^3]
      V         | Volume                                        [m^3]
    \endvartable

    with

    \f[
        E_t = a R_n + b
    \f]

    where
    \vartable
      R_n    | Heat flux through tree crown                    [W/m^2]
      a      | Linear-regression coefficient - slope parameter [g/W/hr]
      b      | Linear-regression coefficient - intercept parameter [g/m^2/hr]
    \endvartable

    Sources applied to either of the below, if exist:
    \verbatim
      T         | Temperature                                [K]
      e         | Internal energy                            [m^2/s^2]
      h         | Enthalphy                                  [m^2/s^2]
    \endverbatim

    Required fields:
    \verbatim
      T         | Temperature                                [K]
      e         | Internal energy                            [m^2/s^2]
      h         | Enthalphy                                  [m^2/s^2]
      LAD       | Leaf area density                          [m^2/m^3]
    \endverbatim

Usage
    Example by using \c constant/fvOptions:
    \verbatim
    evapotranspirationHeatTransfer1
    {
        // Inherited entries
        ...

        // Mandatory entries
        model               tree;

        // Optional entries
        a                   <scalar>;
        b                   <scalar>;
        lambda              <scalar>;
        LAD                 <word>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                         | Type | Reqd | Deflt
      model    | Model name: tree                    | word | yes  | -
      a | Linear-regression coefficient - slope parameter [g/W/hr] | scalar <!--
        -->                                                   | no | 0.3622
      b | Linear-regression coefficient - intercept parameter [g/m^2/hr] <!--
        -->                                         | scalar  | no | 60.758
      lambda   | Specific latent heat of vaporisation [MJ/kg] | scalar <!--
        -->                                                   | no | 2.44
      LAD      | Name of leaf-area density field       | word | no | LAD
    \endtable

    The inherited entries are elaborated in:
      - \link evapotranspirationHeatTransfer.H \endlink
      - \link evapotranspirationHeatTransferModel.H \endlink

SourceFiles
    tree.C
    treeTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_evapotranspirationHeatTransferModels_tree_H
#define Foam_evapotranspirationHeatTransferModels_tree_H

#include "evapotranspirationHeatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace evapotranspirationHeatTransferModels
{

/*---------------------------------------------------------------------------*\
                            Class tree Declaration
\*---------------------------------------------------------------------------*/

class tree
:
    public evapotranspirationHeatTransferModel
{
    // Private Data

        //- Linear-regression coefficient - slope parameter [g/W/hr]
        scalar a_;

        //- Linear-regression coefficient - intercept parameter [g/m^2/hr]
        scalar b_;

        //- Specific latent heat of vaporisation [MJ/kg]
        scalar lambda_;

        //- Name of leaf-area density field
        word LADname_;


    // Private Member Functions

        //- Return area-averaged transpiration [g/m^2/hr]
        //  Calculated from empirical regressions with heat flux
        scalar Et() const;

        //- Return volume weighted leaf area density (LAD) [m^2]
        //  LAD [m^2/m^3]: one-sided area of leaves per unit volume
        tmp<scalarField> leafArea(const labelList& cells) const;


public:

    //- Runtime type information
    TypeName("tree");


    // Constructors

        //- Construct from components
        tree
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        tree(const tree&) = delete;

        //- No copy assignment
        void operator=(const tree&) = delete;


    //- Destructor
    virtual ~tree() = default;


    // Member Functions

    // Evaluation

        //- Return heat-transfer rate for speficied cells [J/s]
        virtual tmp<scalarField> Q(const labelList& cells) const;


    // I-O

        //- Read the dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace evapotranspirationHeatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
