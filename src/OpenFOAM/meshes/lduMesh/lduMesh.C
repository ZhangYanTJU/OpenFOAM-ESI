/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "lduMesh.H"
#include "objectRegistry.H"
#include "processorLduInterface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(lduMesh, 0);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::objectRegistry& Foam::lduMesh::thisDb() const
{
    NotImplemented;
    return NullObjectRef<objectRegistry>();
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<lduMesh>& iproxy
)
{
    const auto& ldum = *iproxy;
    const lduAddressing& addr = ldum.lduAddr();
    const lduInterfacePtrsList interfaces = ldum.interfaces();

    os  << "lduMesh :"
        << " size:" << addr.size()
        << " l:" << addr.lowerAddr().size()
        << " u:" << addr.upperAddr().size()
        << " interfaces:" << interfaces.size()
        << " comm:" << ldum.comm()
        << endl;

    label nCouples = 0;
    forAll(interfaces, i)
    {
        if (interfaces.set(i))
        {
            const labelUList& faceCells = addr.patchAddr(i);
            nCouples += faceCells.size();

            if (isA<processorLduInterface>(interfaces[i]))
            {
                const processorLduInterface& pi = refCast
                <
                    const processorLduInterface
                >(interfaces[i]);

                os  << "    patch:" << i
                    << " type:" << interfaces[i].type()
                    << " size:" << faceCells.size()
                    << " myProcNo:" << pi.myProcNo()
                    << " neighbProcNo:" << pi.neighbProcNo()
                    << " comm:" << pi.comm()
                    << endl;
            }
            else
            {
                os  << "    patch:" << i
                    << " type:" << interfaces[i].type()
                    << " size:" << faceCells.size()
                    << endl;
            }
        }
    }
    os  << "    Interface faces/cells:" << scalar(nCouples)/addr.size()
        << endl;


    // Print actual contents
    if (lduMesh::debug)
    {
        const labelUList& l = addr.lowerAddr();
        const labelUList& u = addr.upperAddr();
        forAll(l, facei)
        {
            os  << "        face:" << facei << " l:" << l[facei]
                << " u:" << u[facei] << endl;
        }
        forAll(interfaces, i)
        {
            if (interfaces.set(i))
            {
                const labelUList& faceCells = addr.patchAddr(i);
                if (faceCells.size())
                {
                    os  << "    patch:" << i
                        << " type:" << interfaces[i].type() << endl;

                    if (isA<processorLduInterface>(interfaces[i]))
                    {
                        const processorLduInterface& pi = refCast
                        <
                            const processorLduInterface
                        >(interfaces[i]);

                        os  << "    myProcNo:" << pi.myProcNo()
                            << " neighbProcNo:" << pi.neighbProcNo()
                            << " comm:" << pi.comm()
                            << endl;
                    }

                    forAll(faceCells, i)
                    {
                        os  << "        " << i << " own:" << faceCells[i]
                            << endl;
                    }
                }
            }
        }
    }

    os.check(FUNCTION_NAME);

    return os;
}


// ************************************************************************* //
