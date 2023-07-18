/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::volVectorField& Foam::objectiveIncompressible::dJdv()
{
    return *dJdvPtr_;
}


inline const Foam::volScalarField& Foam::objectiveIncompressible::dJdp()
{
    return *dJdpPtr_;
}


inline const Foam::volScalarField& Foam::objectiveIncompressible::dJdT()
{
    return *dJdTPtr_;
}


inline const Foam::volScalarField& Foam::objectiveIncompressible::dJdTMvar1()
{
    return *dJdTMvar1Ptr_;
}


inline const Foam::volScalarField& Foam::objectiveIncompressible::dJdTMvar2()
{
    return *dJdTMvar2Ptr_;
}


inline const Foam::fvPatchVectorField&
Foam::objectiveIncompressible::boundarydJdv
(
    const label patchI
)
{
    return bdJdvPtr_()[patchI];
}


inline const Foam::fvPatchScalarField&
Foam::objectiveIncompressible::boundarydJdvn
(
    const label patchI
)
{
    return bdJdvnPtr_()[patchI];
}


inline const Foam::fvPatchVectorField&
Foam::objectiveIncompressible::boundarydJdvt
(
    const label patchI
)
{
    return bdJdvtPtr_()[patchI];
}


inline const Foam::fvPatchVectorField&
Foam::objectiveIncompressible::boundarydJdp
(
    const label patchI
)
{
    return bdJdpPtr_()[patchI];
}


inline const Foam::fvPatchScalarField&
Foam::objectiveIncompressible::boundarydJdT
(
    const label patchI
)
{
    return bdJdTPtr_()[patchI];
}


inline const Foam::fvPatchScalarField&
Foam::objectiveIncompressible::boundarydJdTMvar1
(
    const label patchI
)
{
    return bdJdTMvar1Ptr_()[patchI];
}


inline const Foam::fvPatchScalarField&
Foam::objectiveIncompressible::boundarydJdTMvar2
(
    const label patchI
)
{
    return bdJdTMvar2Ptr_()[patchI];
}


inline const Foam::fvPatchScalarField&
Foam::objectiveIncompressible::boundarydJdnut
(
    const label patchI
)
{
    return bdJdnutPtr_()[patchI];
}


inline const Foam::fvPatchTensorField&
Foam::objectiveIncompressible::boundarydJdGradU
(
    const label patchI
)
{
    return bdJdGradUPtr_()[patchI];
}


inline const Foam::boundaryVectorField&
Foam::objectiveIncompressible::boundarydJdv()
{
    return bdJdvPtr_();
}


inline const Foam::boundaryScalarField&
Foam::objectiveIncompressible::boundarydJdvn()
{
    return bdJdvnPtr_();
}


inline const Foam::boundaryVectorField&
Foam::objectiveIncompressible::boundarydJdvt()
{
    return bdJdvtPtr_();
}


inline const Foam::boundaryVectorField&
Foam::objectiveIncompressible::boundarydJdp()
{
    return bdJdpPtr_();
}


inline const Foam::boundaryScalarField&
Foam::objectiveIncompressible::boundarydJdT()
{
    return bdJdTPtr_();
}


inline const Foam::boundaryScalarField&
Foam::objectiveIncompressible::boundarydJdTMvar1()
{
    return bdJdTMvar1Ptr_();
}


inline const Foam::boundaryScalarField&
Foam::objectiveIncompressible::boundarydJdTMvar2()
{
    return bdJdTMvar2Ptr_();
}


inline const Foam::boundaryScalarField&
Foam::objectiveIncompressible::boundarydJdnut()
{
    return bdJdnutPtr_();
}


inline const Foam::boundaryTensorField&
Foam::objectiveIncompressible::boundarydJdGradU()
{
    return *bdJdGradUPtr_;
}


inline bool Foam::objectiveIncompressible::hasdJdv() const
{
    return bool(dJdvPtr_);
}


inline bool Foam::objectiveIncompressible::hasdJdp() const
{
    return bool(dJdpPtr_);
}


inline bool Foam::objectiveIncompressible::hasdJdT() const
{
    return bool(dJdTPtr_);
}


inline bool Foam::objectiveIncompressible::hasdJdTMVar1() const
{
    return bool(dJdTMvar1Ptr_);
}


inline bool Foam::objectiveIncompressible::hasdJdTMVar2() const
{
    return bool(dJdTMvar2Ptr_);
}


inline bool Foam::objectiveIncompressible::hasBoundarydJdv() const
{
    return bool(bdJdvPtr_);
}


inline bool Foam::objectiveIncompressible::hasBoundarydJdvn() const
{
    return bool(bdJdvnPtr_);
}


inline bool Foam::objectiveIncompressible::hasBoundarydJdvt() const
{
    return bool(bdJdvtPtr_);
}


inline bool Foam::objectiveIncompressible::hasBoundarydJdp() const
{
    return bool(bdJdpPtr_);
}


inline bool Foam::objectiveIncompressible::hasBoundarydJdT() const
{
    return bool(bdJdTPtr_);
}


inline bool Foam::objectiveIncompressible::hasBoundarydJdTMVar1() const
{
    return bool(bdJdTMvar1Ptr_);
}


inline bool Foam::objectiveIncompressible::hasBoundarydJdTMVar2() const
{
    return bool(bdJdTMvar2Ptr_);
}


inline bool Foam::objectiveIncompressible::hasBoundarydJdnut() const
{
    return bool(bdJdnutPtr_);
}


inline bool Foam::objectiveIncompressible::hasBoundarydJdGradU() const
{
    return bool(bdJdGradUPtr_);
}


// ************************************************************************* //
