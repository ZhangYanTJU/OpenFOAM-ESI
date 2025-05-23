/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "forceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

forceList::forceList(liquidFilmBase& film)
:
    PtrList<force>()
{}


forceList::forceList
(
    liquidFilmBase& film,
    const dictionary& dict
)
:
    PtrList<force>()
{
    const wordList models(dict.lookup("forces"));

    Info<< "    Selecting film force models" << endl;
    if (models.size() > 0)
    {
        this->setSize(models.size());

        forAll(models, i)
        {
            set(i, force::New(film, dict, models[i]));
        }
    }
    else
    {
        Info<< "        none" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

forceList::~forceList()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<faVectorMatrix> forceList::correct(areaVectorField& U)
{
    auto tResult = tmp<faVectorMatrix>::New(U, dimForce/dimDensity);
    auto& result = tResult.ref();

    forAll(*this, i)
    {
        result += this->operator[](i).correct(U);
    }

    return tResult;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
