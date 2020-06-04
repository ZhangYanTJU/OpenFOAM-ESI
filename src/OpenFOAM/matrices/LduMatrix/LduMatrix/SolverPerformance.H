/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
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

Class
    Foam::SolverPerformance

Description
    SolverPerformance is the class returned by the LduMatrix solver
    containing performance statistics.

SourceFiles
    SolverPerformance.C

\*---------------------------------------------------------------------------*/

#ifndef SolverPerformance_H
#define SolverPerformance_H

#include "word.H"
#include "FixedList.H"
#include "clockValue.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Type> class SolverPerformance;

template<class Type>
SolverPerformance<Type> max
(
    const SolverPerformance<Type>& sp1,
    const SolverPerformance<Type>& sp2
);

template<class Type>
Istream& operator>>(Istream& is, SolverPerformance<Type>& sp);

template<class Type>
Ostream& operator<<(Ostream& os, const SolverPerformance<Type>& sp);


/*---------------------------------------------------------------------------*\
                       Class SolverPerformance Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class SolverPerformance
{
    // Label type corresponding to Type
    typedef typename pTraits<Type>::labelType labelType;

    // Private Data

        word        solverName_;
        word        fieldName_;
        Type        initialResidual_;
        Type        finalResidual_;
        labelType   nIterations_;
        bool        converged_;
        FixedList<bool, pTraits<Type>::nComponents> singular_;
        FixedList<double, pTraits<Type>::nComponents> timing_;


public:

    // Static Data

        //- Declare name of the class and its debug switch
        ClassName("SolverPerformance");

        //- Large value for the use in solvers
        static constexpr scalar great_ = 1e20;

        //- Small value for the use in solvers
        static constexpr scalar small_ = 1e-20;

        //- Very small value for the use in solvers
        static constexpr scalar vsmall_ = VSMALL;



    // Constructors

        //- Default construct, zero-initialized
        SolverPerformance()
        :
            solverName_(),
            fieldName_(),
            initialResidual_(Zero),
            finalResidual_(Zero),
            nIterations_(Zero),
            converged_(false),
            singular_(false),
            timing_(Zero)
        {}


        //- Construct for given solver and field, zero-initialized
        SolverPerformance
        (
            const word& solverName,
            const word& fieldName,
            const Type& iRes = pTraits<Type>::zero,
            const Type& fRes = pTraits<Type>::zero,
            const labelType& nIter = pTraits<labelType>::zero,
            const bool converged = false,
            const bool singular = false,
            const scalar timingValue = 0
        )
        :
            solverName_(solverName),
            fieldName_(fieldName),
            initialResidual_(iRes),
            finalResidual_(fRes),
            nIterations_(nIter),
            converged_(converged),
            singular_(singular),
            timing_(Zero)
        {
            // Single-value assignment affects single component only
            timing_.first() = timingValue;
        }


    // Member Functions

        //- Return solver name
        const word& solverName() const
        {
            return solverName_;
        }

        //- Return solver name
        word& solverName()
        {
            return solverName_;
        }

        //- Return field name
        const word& fieldName() const
        {
            return fieldName_;
        }


        //- Return initial residual
        const Type& initialResidual() const
        {
            return initialResidual_;
        }

        //- Return initial residual
        Type& initialResidual()
        {
            return initialResidual_;
        }


        //- Return final residual
        const Type& finalResidual() const
        {
            return finalResidual_;
        }

        //- Return final residual
        Type& finalResidual()
        {
            return finalResidual_;
        }


        //- Return number of iterations
        const labelType& nIterations() const
        {
            return nIterations_;
        }

        //- Return number of iterations
        labelType& nIterations()
        {
            return nIterations_;
        }

        //- Has the solver converged?
        bool converged() const
        {
            return converged_;
        }

        //- Is the matrix singular?
        bool singular() const;

        //- Return total timing, the sum of all components
        double timing() const;

        //- Add to single-component timing
        void addTiming(double timingValue)
        {
            timing_.first() += timingValue;
        }

        //- Set single-component timing
        void setTiming(double timingValue)
        {
            timing_.first() = timingValue;
        }

        //- Return the max timing
        double maxTiming() const;

        //- Check, store and return convergence
        bool checkConvergence
        (
            const Type& tolerance,
            const Type& relTolerance
        );

        //- Test singularity
        bool checkSingularity(const Type& residual);

        //- Print summary of solver performance to the given stream
        void print(Ostream& os) const;

        //- Replace component based on the minimal SolverPerformance
        void replace
        (
            const label cmpt,
            const SolverPerformance<typename pTraits<Type>::cmptType>& sp
        );

        //- Return the summary maximum of SolverPerformance<Type>
        //  Effectively it will mostly return solverPerformanceScalar
        SolverPerformance<typename pTraits<Type>::cmptType> max() const;


    // Ostream Operator

        friend Istream& operator>> <Type>
        (
            Istream& is,
            SolverPerformance<Type>& sp
        );

        friend Ostream& operator<< <Type>
        (
            Ostream& os,
            const SolverPerformance<Type>& sp
        );
};


// Global Operators

//- Test for inequality
template<class Type>
bool operator!=
(
    const SolverPerformance<Type>& a,
    const SolverPerformance<Type>& b
);


// Global Functions

//- Return element-wise maximum of SolverPerformance
template<class Type>
SolverPerformance<Type> max
(
    const SolverPerformance<Type>& a,
    const SolverPerformance<Type>& b
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeSolverPerformance(Type)                                           \
    typedef Foam::SolverPerformance<Type> solverPerformance##Type;            \
    defineNamedTemplateTypeNameAndDebug(solverPerformance##Type, 0);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SolverPerformance.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
