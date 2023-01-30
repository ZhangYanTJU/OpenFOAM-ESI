/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
    Foam::multiNodeDecomp

Description
    Decompose given using consecutive application of decomposers,
    to perhaps uneven pieces.
    Note: If uneven pieces are required, the decomposition method
    used must support the processorWeights argument.

SourceFiles
    multiNodeDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef multiNodeDecomp_H
#define multiNodeDecomp_H

#include "decompositionMethod.H"
// #include "List.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class multiNodeDecomp Declaration
\*---------------------------------------------------------------------------*/

class multiNodeDecomp
:
    public decompositionMethod
{
    // Nested classes declarations

    /*---------------------------------------------------------------------------*\
                                Class metaParser Declaration
    \*---------------------------------------------------------------------------*/
    
    // A class responsible for detecting and parsing metadata-related arguments.
    class metaParser {
    public:
        // Detect and return entries related to the given argument.
        // Input:
        //  dict - the coeffs dictionary we are looking inside.
        //  argument - the argument we're looking for.
        //  allowWithoutBrackets - set to true if the argument can be detected without brackets
        //  For example, domains should not be allowed without brackets, but weightsInitialization can.
        static List<string> getEntries(const dictionary& dict, const string& argument, bool allowWithoutBrackets = false);

        // Given a key string of an entry returned from getEntries,
        // Parse and return all ranges described in the key.
        // Note that it is the user's responsibility to make sure that the right endpoint
        // does not exceed the number of children in each node.
        // The user may also write "[]" to specify all children of the node.
        // In this case, the range returned is [0,-1]. Otherwise, endpoints are always non-negative.
        // Input:
        //  key - a key of an entry we are parsing, perhaps from an entry returned in getEntries.
        // Output:
        //  A list of ranges, where the i-th range corresponds to the i-th level of the decomposition tree.
        //  i.e. if two ranges are returned, we will traverse each child of the root in the first range, and recursively
        //  each child in the second range, thus updating nodes at the third level.

        static List<Pair<label>> parseRanges(const string& key);

    };
    enum WeightsInitialization {
            RELATIVE,
            UNIFORM,
            UNKNOWN
        };

    /*---------------------------------------------------------------------------*\
                                Class nodeMetadata Declaration
    \*---------------------------------------------------------------------------*/
    // A class holding all the information necessary for a multi-node decomposition, without building
    // the decompositionMethod objects.
    // The size indicates the number of processors in this subtree. It will be used
    // When computing the offset of the decomposition, and when using the relative weights initialization.
    // The weight is a multiplicative factor applied when decomposing. It is 1 by default and can be set by the user.
    // If the uniform weights initialization is used, all nodes will have the same weight. If the relative weights 
    // initialization is used, each node's weight is set relatively to their size.
    // Then, the weight field can be used to change the weight of a specific node.
    // Note that if the coeffs dictionary contains a processorWeights field, it will not be overwritten.
    // We will then construct a new dictionary with the required numberOfSubdomains and processorWeights.
    
    class nodeMetadata {
    public:
        
        nodeMetadata() : weight(1), size(1), weightsInitialization(UNIFORM), children(0), coeffsDict(nullptr), method(nullptr) {}

        // Constructs a decomposition data tree with dimensions dims and a default method.
        // Input: A list of domains for each level, and a default dictionary method.
        nodeMetadata(const labelList& dims, const dictionary* defaultMethod) : nodeMetadata() {
            initialize(dims, defaultMethod);
        }

        // Initializes an existing nodeMetadata object.
        void initialize(const labelList& dims, const dictionary* defaultMethod) {
            setDict(*defaultMethod);
            constructRecursive(dims, defaultMethod);
        }

        ~nodeMetadata() {
            // Since this class represents a tree, we will need to destruct recursively.
            for(nodeMetadata* child : children)
                delete child;
            
            // Only delete method and dict if they were assigned.
            if(method != nullptr)
                delete method;
            if(coeffsDict != nullptr)
                delete coeffsDict;
        }


        // Getters

        // Get the weight of this node, with respect to the decomposition done in this node's parent
        label getWeight() const {
            return weight;
        }

        // Get the coeffs dictionary for the decomposition of this node.
        const dictionary* getDict() const {
            return coeffsDict;
        }

        // Get the modifiable coeffs dictionary for the decomposition of this node.
        dictionary* getMutableDict() {
            return coeffsDict;
        }

        // Get the number of leaves in this subtree, i.e., number of processors
        // created under this node.
        label getSize() const {
            return size;
        }

        // Get the decomposition method object of this node.
        // Note that construct methods must be called first, otherwise
        // a null pointer will be returned.
        const Foam::decompositionMethod* getMethod() const {
            return method;
        }

        // Get the current weights initialization mode.
        bool getWeightsInitialization() const {
            return weightsInitialization;
        }

        // Get a const pointer to a child of this node.
        const nodeMetadata* getChild(label index) const {
            return children[index];
        }

        // Get a non-const pointer to a child of this node.
        nodeMetadata* getMutableChild(label index) {
            return children[index];
        }

        // Returns the number of direct subdomains this node has.
        label nDomains() const {
            return children.size();
        }
        
        // Returns whether this node represents a leaf (i.e., has no children)
        bool isLeaf() const {
            return children.empty();
        }
        // Setters

        // Set the weight of this node, with respect to the decomposition done in this node's parent
        void setWeight(label weight) {
            this->weight = weight;
        }

        // Set the coeffs dictionary for the decomposition of this node.
        // This creates a copy of the dictionary (and deletes the previous copy created)
        void setDict(const dictionary& dict) {
            if(coeffsDict != nullptr) {
                delete coeffsDict;
            }
            coeffsDict = new dictionary(dict);
        }
        
        // Set the decomposition method object of this node.
        void setMethod(Foam::decompositionMethod* method) {
            this->method = method;
        }

        // Sets the weights initialization mode. If setRecursive is true, propagate to the entire subtree (i.e., the root and all of the descendents)
        void setWeightsInitialization(WeightsInitialization newMode, bool setRecursive = true) {
            weightsInitialization = newMode;
            if(setRecursive) {
                for(nodeMetadata* child : children)
                    child->setWeightsInitialization(newMode, true);
            }
        }

        // Updates
        
        // Update the weights of the nodes at the given indices to the given weight.
        // Input: A string indicating the indices of the nodes to be updated, and
        //        the new weight of the nodes.
        void updateWeight(const string& indices, label newWeight) {
            updateNodes(indices, [newWeight](nodeMetadata* node) {
                node->setWeight(newWeight);
            });
        }
        
        // Update the dimensions array of the nodes at the given indices to the given dimensions array.
        // Input: A string indicating the indices of the nodes to be updated, and
        //        the new list of dimensions.
        void updateDomains(const string& indices,const labelList& dims) {
            updateNodes(indices, [dims](nodeMetadata* node) {
                // Reconstruct using this node's dict.
                // Note that first all domain changes are done,
                // And only then dictionaries are set.
                // So the descendents' dictionaries are not overwritten.
                node->constructRecursive(dims,node->getDict());
            });
        }


        // Update the method of the nodes at the given indices to the given method dictionary.
        // Input: A string indicating the indices of the nodes to be updated, and
        //        the new method dictionary.
        void updateMethod(const string& indices, const dictionary& dict) {
            updateNodes(indices, [dict](nodeMetadata* node) {
                node->setDict(dict);
            });
        }
        
        // Update the weight initialization mode of nodes at the given indices and their descendents to the new mode.
        // Input: A string indicating the indices of the nodes to be updated, and
        //        the new weight mode.
        void updateWeightsInitialization(const string& indices, WeightsInitialization newMode) {
            updateNodes(indices, [newMode](nodeMetadata* node) {
                node->setWeightsInitialization(newMode);
            });
        }

        // Given a list of dictionaries for each level, set the dictionaries accordingly.
        // Input: A list of dictionaries for each level.
        void setLeveledDictionaries(const List<const dictionary*>& dictionaries);

        // To be used within the decompositionMethod's parallelAware function.
        // Returns whether all decompositions in this subtree are parallel aware
        // (i.e., synchronize domains across proc boundaries)
        bool parallelAware() const;

        // Calculate (and add to the dictionary) the new processor weights if reqired,
        // Using the children's weights and the weight initialization mode.
        void updateProcessorWeights();

        // Construct the decompositionMethod object for this node and all its descendents.
        void constructMethods();

    private:

        // The weight of this node in the parent's decomposition, relative to the other nodes.
        // Overrides weights set by the weights initialization.
        label weight;

        // The size of the node indicates the total number of subdomains this node has.
        label size;
        
        // An enum describing the weights initialization.
        WeightsInitialization weightsInitialization;
        
        // The direct descendents.
        List<nodeMetadata*> children;
        
        // The dictionary used to construct the decomposition method.
        dictionary* coeffsDict; 
        
        // The decomposition method of this node.
        const Foam::decompositionMethod* method;

        // Recursively constructs the subtree rooted at this node
        // Input: A list of dimensions and the dictionary of the default method.
        void constructRecursive(const labelList& dims, const dictionary* defaultMethod);

        // Update all nodes matching the given indices with the given updating function.
        // Input: A list of ranges for each level, and a function that receives a pointer to nodeMetadata, 
        // updates it accordingly and returns nothing.
        void updateNodes(const string& key, const std::function<void(nodeMetadata*)>& update);

        // Internal implementation of updateNodes.
        // The list of ranges are constructed by passing the key argument to the meta parser.
        // nCurIdx is an internal variable that indicates our location inside the indices array.
        void updateNodes(const List<Pair<label>>& indices, const std::function<void(nodeMetadata*)>& update, label nCurIdx = 0);

        // This function is used inside the public setLeveledDictionaries function.
        void setLeveledDictionaries(const List<const dictionary*>& dictionaries, label nLevel);

        // Parse the range of indices starting at the (string) index nStartIndex.
        // Input: The indices string, and the starting position of the range
        // (i.e, the position of the opening bracket)
        // Returns a pair representing the range if succeeded,
        // or crashes the program with an appropriate error message if failed to parse.
        Pair<label> parseRange(const string& indices, label nStartIndex) const;

    };


    // Private Data

        //- The decomposition metadata.
        nodeMetadata rootMetadata_;


    // Private Member Functions

        //- Read coeffsDict and construct the decomposition metadata.
        void initializeMetadata(const dictionary& coeffsDict);


        //- Given connectivity across processors work out connectivity
        //  for a (consistent) subset
        void subsetGlobalCellCells
        (
            const label nDomains,
            const label domainI,
            const labelList& dist,

            const labelListList& cellCells,
            const labelList& set,
            labelListList& subCellCells,
            labelList& cutConnections
        ) const;

        //- Decompose at 'currLevel' without addressing
        void decompose
        (
            const labelListList& pointPoints,
            const pointField& points,
            const scalarField& pointWeights,
            const labelUList& pointMap,  // map back to original points
            const nodeMetadata& decomposeData,
            const label leafOffset,

            labelList& finalDecomp
        ) const;


        //- No copy construct
        multiNodeDecomp(const multiNodeDecomp&) = delete;

        //- No copy assignment
        void operator=(const multiNodeDecomp&) = delete;

public:

    //- Runtime type information
    TypeName("multiNode");


    // Constructors

        //- Construct given decomposition dictionary and optional region name
        explicit multiNodeDecomp
        (
            const dictionary& decompDict,
            const word& regionName = ""
        );


    //- Destructor
    virtual ~multiNodeDecomp() = default;


    // Member Functions

        //- Is method parallel aware?
        //  i.e. does it synchronize domains across proc boundaries
        virtual bool parallelAware() const;

        //- Inherit decompose from decompositionMethod
        using decompositionMethod::decompose;

        //- Return for every coordinate the wanted processor number.
        //  Use the mesh connectivity (if needed)
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& points,
            const scalarField& pointWeights
        ) const;

        //- Return for every coordinate the wanted processor number.
        //  Explicitly provided connectivity - does not use mesh_.
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cc,
            const scalarField& cWeights
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
