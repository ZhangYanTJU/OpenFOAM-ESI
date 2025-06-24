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
    Foam::functionObjects::fieldStatistics

Group
    grpFieldFunctionObjects

Description
    Calculates various statistics of the specified fields.

    Operands:
    \table
      Operand      | Type                          | Location
      input        | vol\<Type\>Field(s)           | Registry
      output file  | dat | \<case\>/postProcessing/\<FO\>/\<time\>/\<file\>
      output field | -                             | -
    \endtable

    where \c \<Type\>=Scalar/Vector/SphericalTensor/SymmTensor/Tensor.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    fieldStatistics1
    {
        // Mandatory entries
        type             fieldStatistics;
        libs             (fieldFunctionObjects);
        fields           (<wordList>);
        statistics       (<wordList>);

        // Optional entries
        mode             <word>;
        mean             <word>;
        extrema          <bool>;
        internal         <bool>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd  | Deflt
      type         | Type name: fieldStatistics         | word | yes   | -
      libs         | Library name: fieldFunctionObjects | word | yes   | -
      fields       | List of operand fields             | wordList | yes | -
      statistics   | List of operand statistics         | wordList | yes | -
      mode   | Output format of the statistical results | word | no | magnitude
      mean     | Type of the mean operation             | word | no | arithmetic
      internal | Flag to use internal fields only in computing statistics <!--
               -->                                      | bool | no | false
      extrema | Flag to enable extrema data calculations | bool | no | false
    \endtable

    Available statistics of the operand field through the \c statistics entry:
    \verbatim
      min        | Minimum value
      max        | Maximum value
      mean       | Arithmetic mean value (optionally volume-weighted)
      variance   | Sample variance value (unbiased)
    \endverbatim

    Options for the \c mode entry:
    \verbatim
      magnitude  | Output statistics magnitude-wise
      component  | Output statistics separately for each component
    \endverbatim

    Options for the \c mean entry:
    \verbatim
      arithmetic  | Arithmetic mean (average)
      volumetric  | Volume-weighted arithmetic mean
    \endverbatim

    The inherited entries are elaborated in:
      - \link functionObject.H \endlink
      - \link writeFile.H \endlink

SourceFiles
    fieldStatistics.cxx
    fieldStatisticsImpl.cxx

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fieldStatistics_H
#define functionObjects_fieldStatistics_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldSelection.H"
#include <functional>
#include <variant>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class fieldStatistics Declaration
\*---------------------------------------------------------------------------*/

class fieldStatistics
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private Enumerations

        //- Options for the output format of the statistical results
        enum class modeType : char
        {
            MAG = 0,    //!< "Output statistics magnitude-wise"
            CMPT        //!< "Output statistics separately for each component"
        };

        //- Names for modeType
        static const Enum<modeType> modeTypeNames_;

        //- Options for the mean type of the specified fields
        enum class meanType : char
        {
            ARITHMETIC = 0,  //!< "Arithmetic mean (average)"
            VOLUMETRIC       //!< "Volume-weighted arithmetic mean"
        };

        //- Names for meanType
        static const Enum<meanType> meanTypeNames_;

        //- Options for the type of statistics calculation
        enum class calcType : char
        {
            UNKNOWN = 0,    //!< placeholder
            MIN,            //!< min value
            MAX,            //!< max value
            MEAN,           //!< mean value
            VARIANCE        //!< variance
        };

        //- Names for calcType
        static const Enum<calcType> calcTypeNames_;


    // Private Classes

        //- Type-safe union for input field types
        using variantInput = std::variant
        <
            scalarField,
            vectorField,
            sphericalTensorField,
            symmTensorField,
            tensorField
        >;

        //- Type-safe union for output data types
        using variantOutput = std::variant
        <
            scalar,
            vector,
            sphericalTensor,
            symmTensor,
            tensor
        >;

        //- Class to encapsulate information about specified statistic
        struct statistic
        {
            //- Name of the statistic
            word name_;

            //- Returns the value of the specified statistic
            std::function<variantOutput(variantInput)> calc;
        };

        //- Class to encapsulate the data about minima and maxima
        struct extremaData
        {
            //- Value of the extremum
            variantOutput value_;

            //- Processor index of the extremum
            label procID_;

            //- Cell index of the extremum
            label cellID_;

            //- Position (cell or face centre) of the extremum
            point position_;
        };


    // Private Data

        //- Flag to use internal fields only in computing statistics
        bool internal_;

        //- Flag to enable extrema data calculations
        bool extrema_;

        //- Output-format mode - only applicable for tensor ranks > 0
        modeType mode_;

        //- Type of the mean of the specified fields
        meanType mean_;

        //- Operand fields on which statistics are computed
        volFieldSelection fieldSet_;

        //- List of file pointers; one file per field
        HashPtrTable<OFstream> filePtrs_;

        //- List of file pointers for extrema data; one file per field
        HashPtrTable<OFstream> extremaFilePtrs_;

        //- Hash table containing all specified statistics
        HashTable<statistic> statistics_;

        //- Hash table containing all statistical results per field
        HashTable<HashTable<variantOutput>> results_;

        //- Hash table containing the results of the extrema per field
        HashTable<Pair<extremaData>> extremaResults_;


    // Private Member Functions

        //- Return the statistic container
        statistic createStatistic(const word& statName, const modeType mode);

        //- Compute the specified statistics of a given field
        template<class T>
        bool calcStat(const word& fieldName);


        // Central tendency statistics

        //- Return the arithmetic mean of the given input field
        template<class T>
        T calcMean(const Field<T>& field) const;


        // Dispersion statistics

        //- Return the minimum value of the given input field
        //  Store the processor index, cell index and location of the minimum
        template<class T>
        T calcMin(const Field<T>& field) const;

        //- Return the maximum value of the given input field
        //  Store the processor index, cell index and location of the maximum
        template<class T>
        T calcMax(const Field<T>& field) const;

        //- Return the sample variance of the given input field
        template<class T>
        T calcVariance(const Field<T>& field) const;

        //- Return a combined field: internal + flattened boundary
        template<class GeoField>
        tmp<Field<typename GeoField::value_type>>
        flatten(const GeoField& field) const;

        //- Return the extrema data of the specified field
        template<class GeoField>
        Pair<extremaData> calcExtremaData(const GeoField& field) const;

        //- Output the file header information
        void writeFileHeader(Ostream& os, const word& fieldName);

        //- Write the statistical data to the specified file
        void writeStatData();

        //- Write the statistical data to the standard stream
        void logStatData();

        //- Output the extrema-data file header information
        void writeExtremaFileHeader(Ostream& os, const word& fieldName);

        //- Write extrema data to the specified file
        void writeExtremaData();

        //- Write extrema data to the standard stream
        void logExtremaData();


public:

    //- Runtime type information
    TypeName("fieldStatistics");


    // Generated Methods

        //- No copy construct
        fieldStatistics(const fieldStatistics&) = delete;

        //- No copy assignment
        void operator=(const fieldStatistics&) = delete;


    // Constructors

        //- Construct from name, Time and dictionary
        fieldStatistics
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~fieldStatistics() = default;


    // Member Functions

        //- Read function object settings
        virtual bool read(const dictionary&);

        //- Execute function object
        virtual bool execute();

        //- Write function object results
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
