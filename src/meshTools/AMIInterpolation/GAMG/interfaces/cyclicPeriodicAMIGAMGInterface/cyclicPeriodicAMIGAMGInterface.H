/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::cyclicPeriodicAMIGAMGInterface

Description
    GAMG agglomerated cyclic periodic AMI interface.

SourceFiles
    cyclicPeriodicAMIGAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicPeriodicAMIGAMGInterface_H
#define cyclicPeriodicAMIGAMGInterface_H

#include "GAMGInterface.H"
#include "cyclicPeriodicAMILduInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class cyclicPeriodicAMIGAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class cyclicPeriodicAMIGAMGInterface
:
    public GAMGInterface,
    virtual public cyclicPeriodicAMILduInterface
{
    // Private data

        //- Reference for the cyclicLduInterface from which this is
        //  agglomerated
        const cyclicPeriodicAMILduInterface& fineCyclicPeriodicAMIInterface_;

        //- AMI interface
        autoPtr<AMIPatchToPatchInterpolation> amiPtr_;


    // Private Member Functions

        //- No copy construct
        cyclicPeriodicAMIGAMGInterface(const cyclicPeriodicAMIGAMGInterface&) =
            delete;

        //- No copy assignment
        void operator=(const cyclicPeriodicAMIGAMGInterface&) = delete;


public:

    //- Runtime type information
    TypeName("cyclicPeriodicAMI");


    // Constructors

        //- Construct from fine level interface,
        //- local and neighbour restrict addressing
        cyclicPeriodicAMIGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const lduInterface& fineInterface,
            const labelField& restrictAddressing,
            const labelField& neighbourRestrictAddressing,
            const label fineLevelIndex,
            const label coarseComm
        );


    //- Destructor
    virtual ~cyclicPeriodicAMIGAMGInterface();


    // Member Functions

        // Interface transfer functions

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const;


        //- Cyclic interface functions

            //- Return neighbour processor number
            virtual label neighbPatchID() const
            {
                return fineCyclicPeriodicAMIInterface_.neighbPatchID();
            }

            virtual bool owner() const
            {
                return fineCyclicPeriodicAMIInterface_.owner();
            }

            virtual const cyclicPeriodicAMIGAMGInterface& neighbPatch() const
            {
                return dynamic_cast<const cyclicPeriodicAMIGAMGInterface&>
                (
                    coarseInterfaces_[neighbPatchID()]
                );
            }

            virtual const AMIPatchToPatchInterpolation& AMI() const
            {
                return *amiPtr_;
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return fineCyclicPeriodicAMIInterface_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return fineCyclicPeriodicAMIInterface_.reverseT();
            }


        // I/O

            //- Write to stream
            virtual void write(Ostream&) const
            {
                //TBD. How to serialise the AMI such that we can stream
                // cyclicPeriodicAMIGAMGInterface.
                NotImplemented;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
