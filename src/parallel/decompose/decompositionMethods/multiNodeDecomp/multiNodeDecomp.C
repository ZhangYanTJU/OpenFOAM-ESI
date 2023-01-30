/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "multiNodeDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "globalIndex.H"
#include "mapDistribute.H"
#include "DynamicList.H"





// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiNodeDecomp, 0);
    addToRunTimeSelectionTable
    (
        decompositionMethod,
        multiNodeDecomp,
        dictionary
    );
}



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
namespace Foam {
    void multiNodeDecomp::initializeMetadata(const dictionary& coeffsDict) {
        word defaultMethod;
        dictionary defaultMethodDict;
        if(coeffsDict.readIfPresent("method", defaultMethod, keyType::LITERAL)) {
            defaultMethodDict.add("method",defaultMethod);
            const dictionary& subMethodCoeffsDict
            (
                findCoeffsDict
                (
                    coeffsDict,
                    defaultMethod + "Coeffs",
                    selectionType::NULL_DICT
                )
            );
            if(subMethodCoeffsDict.size())
                defaultMethodDict.add(subMethodCoeffsDict.dictName(), subMethodCoeffsDict);
        }
        labelList domains;

        label nTotal = 0;
        label nLevels = 0;

        //Check if any meta argument is changed using the new syntax.
        //If they are, we cannot infer an additional level of decomposition,
        //as it may interfere with the indices.
        List<string> domainChanges = metaParser::getEntries(coeffsDict, "domains");
        List<string> methodChanges = metaParser::getEntries(coeffsDict, "method");
        List<string> weightChanges = metaParser::getEntries(coeffsDict, "weight");
        //We can parse weightMode without brackets too
        List<string> weightModeChanges = metaParser::getEntries(coeffsDict, "weightsInitialization", true);

        bool bChangesDomains = !domainChanges.empty();
        bool bChangesArguments =        bChangesDomains
                                    || (!methodChanges.empty())
                                    || (!weightChanges.empty())
                                    || (!weightModeChanges.empty());
            
        bool bMetadataInitialized = false;

        // Found (non-recursive, no patterns) "method" and "domains" ?
        // Allow as quick short-cut entry
        if
        (
            // non-recursive, no patterns
            coeffsDict.readIfPresent("method", defaultMethod, keyType::LITERAL)
            // non-recursive, no patterns
        && coeffsDict.readIfPresent("domains", domains, keyType::LITERAL)
        )
        {
            // Short-cut version specified by method, domains only

            nTotal = (domains.empty() ? 0 : 1);
            for (const label n : domains)
            {
                nTotal *= n;
                ++nLevels;
            }
            
            //Update domains here
            if(nTotal != 0 && bChangesDomains) {
                rootMetadata_.initialize(
                    domains,
                    &defaultMethodDict
                );
                bMetadataInitialized = true;
                for(string key : domainChanges)
                    rootMetadata_.updateDomains(   key,
                                                    coeffsDict.get<labelList>(key, keyType::LITERAL));
                
                nTotal = rootMetadata_.getSize();
            }

            if (nTotal == 1)
            {
                // Emit Warning
                nTotal = nDomains();
                nLevels = 1;
                domains.setSize(1);
                domains[0] = nTotal;
            }
            //If bChangesDomains is true, we do not want to add another dimension as this
            //may affect the user's assignments of domains/weights/methods later on.
            else if (nTotal > 0 && nTotal < nDomains() && !(nDomains() % nTotal) && !bChangesArguments)
            {
                // nTotal < nDomains, but with an integral factor,
                // which we insert as level 0
                ++nLevels;

                labelList old(std::move(domains));

                domains.setSize(old.size()+1);

                domains[0] = nDomains() / nTotal;
                forAll(old, i)
                {
                    domains[i+1] = old[i];
                }
                nTotal *= domains[0];


                Info<<"    inferred level 0 with " << domains[0]
                    << " domains" << nl << nl;
            }

            if (!nLevels || nTotal != nDomains())
            {
                FatalErrorInFunction
                    << "Top level decomposition specifies " << nDomains()
                    << " domains which is not equal to the product of"
                    << " all sub domains " << nTotal
                    << exit(FatalError);
            }

            if(!bMetadataInitialized) {
                bMetadataInitialized = true;
                rootMetadata_.initialize(
                    domains,
                    &defaultMethodDict
                );
            }
        }
        else
        {
            // Specified by full dictionaries

            // Create editable methods dictionaries
            // - Only consider sub-dictionaries with a "numberOfSubdomains" entry
            //   This automatically filters out any coeffs dictionaries

            label nTotal = 1;
            List<const dictionary*> methods;
            for (const entry& dEntry : coeffsDict)
            {
                word methodName;

                if
                (
                    dEntry.isDict()
                    // non-recursive, no patterns
                && dEntry.dict().found("numberOfSubdomains", keyType::LITERAL)
                )
                {
                    domains.append(dEntry.dict().get<label>("numberOfSubdomains"));
                    nTotal *= domains.last();
                    // No method specified? can use a default method?

                    const bool addDefaultMethod
                    (
                        !(dEntry.dict().found("method", keyType::LITERAL))
                    && !defaultMethod.empty()
                    );
                    if(!(dEntry.dict().found("method",keyType::LITERAL)) && defaultMethod.empty()) {
                        FatalErrorInFunction <<
                            dEntry.keyword() <<
                            " dictionary does not contain method, and no default method is specified."
                            << nl << exit(FatalError);
                    }
                    dictionary* levelDict = new dictionary(dEntry.dict());
                    levelDict->remove("numberOfSubdomains");
                    if(addDefaultMethod) levelDict->add("method", defaultMethod);
                    methods.append(levelDict);
                }
            }
            if(domains.empty())
                nTotal = 0;


            rootMetadata_.initialize(domains, methods[0]);
            bMetadataInitialized = true;
            for(string key : domainChanges)
                rootMetadata_.updateDomains(   key,
                                                coeffsDict.get<labelList>(key, keyType::LITERAL));
            
            if(nTotal != nDomains()) {
                FatalErrorInFunction
                    << "Top level decomposition specifies " << nDomains()
                    << " domains which is not equal to the product of"
                    << " all sub domains " << nTotal << " manually defined by dictionaries. "
                    << exit(FatalError);
            }
            rootMetadata_.setLeveledDictionaries(methods);
            for(const dictionary* method : methods)
                delete method;
        }


        for(string key : methodChanges)
            rootMetadata_.updateMethod(key, coeffsDict.subDict(key, keyType::LITERAL));
        
        for(string key : weightChanges)
            rootMetadata_.updateWeight(key, coeffsDict.get<label>(key, keyType::LITERAL));

        for(string key : weightModeChanges) {
            word value = coeffsDict.get<word>(key, keyType::LITERAL);
            WeightsInitialization newValue = UNKNOWN;
            
            if(value=="uniform")
                newValue = UNIFORM;
            else if(value == "relative")
                newValue = RELATIVE;
            else
                FatalErrorInFunction <<
                    "unknown weights initialization (" << value << "). Must be one of: relative, uniform."
                    << nl << exit(FatalError);
            
            rootMetadata_.updateWeightsInitialization(key, newValue);
        }
        
        if(!rootMetadata_.isLeaf())
            rootMetadata_.constructMethods();
    }


    // Given a subset of cells determine the new global indices. The problem
    // is in the cells from neighbouring processors which need to be renumbered.
    void multiNodeDecomp::subsetGlobalCellCells
    (
        const label nDomains,
        const label domainI,
        const labelList& dist,

        const labelListList& cellCells,
        const labelList& set,
        labelListList& subCellCells,
        labelList& cutConnections
    ) const
    {
        // Determine new index for cells by inverting subset
        labelList oldToNew(invert(cellCells.size(), set));

        globalIndex globalCells(cellCells.size());

        // Subset locally the elements for which I have data
        subCellCells = UIndirectList<labelList>(cellCells, set);

        // Get new indices for neighbouring processors
        List<Map<label>> compactMap;
        mapDistribute map(globalCells, subCellCells, compactMap);
        map.distribute(oldToNew);
        labelList allDist(dist);
        map.distribute(allDist);

        // Now we have:
        // oldToNew : the locally-compact numbering of all our cellCells. -1 if
        //            cellCell is not in set.
        // allDist  : destination domain for all our cellCells
        // subCellCells : indexes into oldToNew and allDist

        // Globally compact numbering for cells in set.
        globalIndex globalSubCells(set.size());

        // Now subCellCells contains indices into oldToNew which are the
        // new locations of the neighbouring cells.

        cutConnections.setSize(nDomains);
        cutConnections = 0;

        forAll(subCellCells, subCelli)
        {
            labelList& cCells = subCellCells[subCelli];

            // Keep the connections to valid mapped cells
            label newI = 0;
            forAll(cCells, i)
            {
                // Get locally-compact cell index of neighbouring cell
                const label nbrCelli = oldToNew[cCells[i]];
                if (nbrCelli == -1)
                {
                    cutConnections[allDist[cCells[i]]]++;
                }
                else
                {
                    // Reconvert local cell index into global one

                    // Get original neighbour
                    const label celli = set[subCelli];
                    const label oldNbrCelli = cellCells[celli][i];
                    // Get processor from original neighbour
                    const label proci = globalCells.whichProcID(oldNbrCelli);
                    // Convert into global compact numbering
                    cCells[newI++] = globalSubCells.toGlobal(proci, nbrCelli);
                }
            }
            cCells.setSize(newI);
        }
    }


    void multiNodeDecomp::decompose
    (
        const labelListList& pointPoints,
        const pointField& points,
        const scalarField& pointWeights,
        const labelUList& pointMap,     // map back to original points
        const nodeMetadata& decomposeData,
        const label leafOffset,

        labelList& finalDecomp
    ) const
    {
        labelList dist
        (
            decomposeData.getMethod()->decompose
            (
                pointPoints,
                points,
                pointWeights
            )
        );

        // Number of domains at the current level
        const label nCurrDomains = decomposeData.nDomains();

        // Calculate the domain remapping.
        // The decompose() method delivers a distribution of [0..nDomains-1]
        // which we map to the final location according to the decomposition
        // leaf we are on.

        labelList domainOffsets(nCurrDomains);
        domainOffsets[0] = leafOffset;
        for(label nDomain = 1; nDomain < nCurrDomains; ++nDomain) {
            domainOffsets[nDomain] = domainOffsets[nDomain-1] + decomposeData.getChild(nDomain-1)->getSize();
        }

        // Extract processor+local index from point-point addressing
        forAll(pointMap, i)
        {
            finalDecomp[pointMap[i]] = domainOffsets[dist[i]];
        }

        if (nCurrDomains > 0)
        {
            // Recurse

            // Determine points per domain
            labelListList domainToPoints(invertOneToMany(nCurrDomains, dist));

            for (label domainI = 0; domainI < nCurrDomains; ++domainI)
            {
                if(decomposeData.getChild(domainI)->isLeaf()) continue;
                // Extract elements for current domain
                const labelList domainPoints(findIndices(dist, domainI));

                // Subset point-wise data.
                pointField subPoints(points, domainPoints);
                scalarField subWeights(pointWeights, domainPoints);
                labelList subPointMap(labelUIndList(pointMap, domainPoints));
                // Subset point-point addressing (adapt global numbering)
                labelListList subPointPoints;
                labelList nOutsideConnections;
                subsetGlobalCellCells
                (
                    nCurrDomains,
                    domainI,
                    dist,

                    pointPoints,
                    domainPoints,

                    subPointPoints,
                    nOutsideConnections
                );

                decompose
                (
                    subPointPoints,
                    subPoints,
                    subWeights,
                    subPointMap,
                    *decomposeData.getChild(domainI),
                    domainOffsets[domainI], // The offset for this level and leaf

                    finalDecomp
                );
            }
        }
    }


    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    multiNodeDecomp::multiNodeDecomp
    (
        const dictionary& decompDict,
        const word& regionName
    )
    :
        decompositionMethod(decompDict, regionName),
        rootMetadata_()
    {
        const dictionary& coeffsDict(
            findCoeffsDict(
                typeName + "Coeffs",
                (selectionType::EXACT | selectionType::MANDATORY)
            )
        );
        initializeMetadata(coeffsDict);
    }

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


    bool multiNodeDecomp::parallelAware() const
    {
        return rootMetadata_.parallelAware();
    }


    labelList multiNodeDecomp::decompose
    (
        const polyMesh& mesh,
        const pointField& cc,
        const scalarField& cWeights
    ) const
    {
        CompactListList<label> cellCells;
        calcCellCells(mesh, identity(cc.size()), cc.size(), true, cellCells);

        labelList finalDecomp(cc.size(), Zero);
        labelList cellMap(identity(cc.size()));

        decompose
        (
            cellCells.unpack(),
            cc,
            cWeights,
            cellMap,      // map back to original cells
            rootMetadata_,
            0,

            finalDecomp
        );

        return finalDecomp;
    }


    labelList multiNodeDecomp::decompose
    (
        const labelListList& globalPointPoints,
        const pointField& points,
        const scalarField& pointWeights
    ) const
    {
        labelList finalDecomp(points.size(), Zero);
        labelList pointMap(identity(points.size()));

        decompose
        (
            globalPointPoints,
            points,
            pointWeights,
            pointMap,       // map back to original points
            rootMetadata_,
            0,

            finalDecomp
        );

        return finalDecomp;
    }

// * * * * * * * * * * * * Meta Parser Class * * * * * * * * * * * * * //



    List<string> multiNodeDecomp::metaParser::getEntries(const dictionary& dict, const string& argument, bool allowWithoutBrackets) {
        string argumentBracket = argument + "[";
        DynamicList<string, 4> Result;
        for(auto& dEntry : dict) {
            if(dEntry.keyword().starts_with(argumentBracket) || (allowWithoutBrackets && dEntry.keyword() == argument))
                Result.push_back(dEntry.keyword());
        }
        return Result;
    }

    List<Pair<label>> multiNodeDecomp::metaParser::parseRanges(const string& key) {
        
        // First, discard the argument and process the indices only.
        // The current syntax is argument[...]
        // Assuming that this key was returned in getEntries,
        // if there is no '[', it is OK and we use the
        // empty string (update the root). 
        string indices = "";
        if(key.find_first_of('[') != key.npos) {
            // There is a '[' in the string.
            // We can substr from that location.
            label nFirstBracket = key.find('[');    
            indices = key.substr(nFirstBracket);
        }

        // All checks print an error message if failed, explaining why.

        DynamicList<Pair<label>, 4> Result;
        label nCurPtr = 0, nIndicesLength = indices.size();
        // As long as there are more ranges to parse.
        while(nCurPtr != nIndicesLength) {
            // First, check if there is an opening bracket.
            if(indices[nCurPtr]!='[')
                FatalError
                    << "Error when parsing indices "
                    << indices << ": Expected '[', found "
                    << indices[nCurPtr] << ". Aborting\n"
                    << exit(FatalError);
            
            // Then, find the matching close bracket.
            label nEndIndex = indices.find(']', nCurPtr);
            if(nEndIndex == nIndicesLength) {
                FatalError
                    << "Error when parsing indices "
                    << indices << ": Expected ']' after '['. Aborting\n"
                    << exit(FatalError);
            }
            // Read inside the brackets, mark the hyphen if it exists, and make sure
            // every character is either a digit or a hyphen.
            // Note that only one hyphen may exist.
            label nHyphenIdx=-1;
            for(label nCurIndex = nCurPtr+1; nCurIndex < nEndIndex; ++nCurIndex) {
                if(!isdigit(indices[nCurIndex])&&indices[nCurIndex]!='-') {
                        FatalError
                            << "Error when parsing indices "
                            << indices << ": Expected digit/'-'/']', found "
                            << indices[nCurIndex] << ". Aborting\n"
                            << exit(FatalError);
                }
                if(indices[nCurIndex]=='-') {
                    if(nHyphenIdx!=-1)
                        FatalError
                        << "Error when parsing indices "
                        << indices << ": Found two hyphens(-) inside an index. Aborting\n"
                        << exit(FatalError);

                    nHyphenIdx = nCurIndex;
                }
            }
            label nLeft,nRight;
            if(nHyphenIdx == -1) {
                // Not a range - just a single index, or empty brackets (indicating to change the whole range).
                if(nCurPtr+1==nEndIndex) nLeft = 0, nRight = -1;
                else {
                    string sNum = indices.substr(nCurPtr+1,nEndIndex-nCurPtr-1);
                    nLeft = nRight = atoi(sNum.c_str());
                }
            } else {
                // A range of indices.
                // Assert that the hyphen is not right next to the brackets.
                if(nHyphenIdx+1==nEndIndex||nCurPtr+1==nHyphenIdx)
                    FatalError
                        << "Error when parsing indices "
                        << indices << ": Expected number, found "
                        << (nCurPtr+1==nHyphenIdx?'-':']')
                        << ". Aborting\n"
                        << exit(FatalError);

                // Parse the numbers
                string sLeftNum = indices.substr(nCurPtr+1,nHyphenIdx-nCurPtr-1);
                string sRightNum = indices.substr(nHyphenIdx+1,nEndIndex-nHyphenIdx-1);
                nLeft = atoi(sLeftNum.c_str());
                nRight = atoi(sRightNum.c_str());
                // Make sure left endpoint is at most the right endpoint
                if(nLeft>nRight)
                    FatalError
                    << "Error when parsing indices "
                    << indices << ": right endpoint("<< nRight
                    << ") cannot be smaller than left endpoint("
                    << nLeft << "). Aborting\n"
                    << exit(FatalError);
            }
            // Move the pointer after the closing bracket and append to the result list.
            nCurPtr = nEndIndex + 1;
            Result.push_back({nLeft,nRight});
        }
        return Result;        
    }
    
// * * * * * * * * * * * * Node Metadata Class * * * * * * * * * * * * * //

    void multiNodeDecomp::nodeMetadata::setLeveledDictionaries(const List<const dictionary*>& dictionaries) {
        setLeveledDictionaries(dictionaries, 0);
    }

    bool multiNodeDecomp::nodeMetadata::parallelAware() const {
        // The decomposition tree is parallel aware if and only if all methods used are parallel aware.
        // If this is a leaf, we are OK.
        if(children.empty())
            return true;

        // Otherwise, check if the method used in this node is parallel aware.
        if(!method->parallelAware())
            return false;
        
        // Check recursively, and if any child is not parallel aware - return false.
        for(auto& child : children)
            if(!child->parallelAware())
                return false;
        
        return true;
    }


    void multiNodeDecomp::nodeMetadata::updateProcessorWeights() {
        label nDom = nDomains();
        word methodCoeffsName = coeffsDict->get<word>("method") + "Coeffs";
        // If processorWeights were set by the user, we do not modify them.
        if(
            // Check if the user did not specify processorWeights under the coeffs dictionary or the methodCoeffs dictionary
            !(coeffsDict->subDictOrAdd(methodCoeffsName).found("processorWeights", keyType::LITERAL)
            || coeffsDict->subDictOrAdd("coeffs").found("processorWeights", keyType::LITERAL))) {
            // Then we should compute weights on our own
            Field<float> processorWeights(nDom);
            forAll(children, i) {
                if(children[i]->weight != 1)
                    processorWeights[i] = children[i]->weight;
                else switch(weightsInitialization) {
                    case RELATIVE:
                        processorWeights[i] = children[i]->size;
                        break;
                    case UNIFORM:
                        processorWeights[i] = 1;
                        break;
                    default:
                        FatalError
                            << "Weights initialization is not handled in updateProcessorWeights. Aborting\n"
                            << exit(FatalError);
                }
            }
            
            coeffsDict->subDictOrAdd(methodCoeffsName).add("processorWeights", processorWeights);
        }
    }
    void multiNodeDecomp::nodeMetadata::constructMethods() {
        // Special handling of nDomains = 1, because some decomposition methods crash when decomposing to one domain.
        label nDom = nDomains();
        if(nDom==1) {
            coeffsDict->clear();
            coeffsDict->add("method","none");
        } else
            updateProcessorWeights();
        coeffsDict->add("numberOfSubdomains",nDom);

        // Non-verbose construction of decomposition methods would be nice
        method = decompositionMethod::New(*coeffsDict).release();
        // Cannot release coeffsDict from memory because method uses a reference that must stay alive

        forAll(children, i) {
            if(!children[i]->isLeaf())
                children[i]->constructMethods();
        }
    }

    // Recursively construct the decomposition tree, given the list of dimensions and a default method.
    void multiNodeDecomp::nodeMetadata::constructRecursive(const labelList& dims, const dictionary* defaultMethod) {
        if(!dims.empty()) {
            // The list of dimensions of the children is the current list without the first element.
            labelList newDims(dims.size() - 1);
            forAll(newDims, i)
                newDims[i] = dims[i+1];
            
            // Construct children recursively
            // First, resize existing children
            // And delete the excess
            forAll(children, i) {
                if(i < dims[0]) 
                    children[i]->constructRecursive(newDims, defaultMethod);
                else
                    delete children[i];
            }
            label nOldSize = children.size();
            children.resize(dims[0]);
            // If the new array is bigger we will need to allocate new children.
            for(label i = nOldSize; i < dims[0]; ++i)
                children[i] = new nodeMetadata(newDims, defaultMethod);
            
            // Compute size (number of leaves in subtree)
            size = dims[0];
            if(!children.empty())
                size *= children[0]->size;
        }
    }
    void multiNodeDecomp::nodeMetadata::updateNodes(const string& key, const std::function<void(nodeMetadata*)>& update) {
        List<Pair<label>> indicesList = metaParser::parseRanges(key);
        updateNodes(indicesList, update);
    }

    // Parse the indices, and apply the update function to all matching nodes.
    // nCurPtr is used to indicate the index we are now parsing (instead of sending substrings of indices)
    void multiNodeDecomp::nodeMetadata::updateNodes(const List<Pair<label>>& indices, const std::function<void(nodeMetadata*)>& update, label nCurIdx) {
        if(nCurIdx == label(indices.size())) update(this);
        else {
            // Otherwise, call recursively.
            label nLeft, nRight, nChildren = children.size();
            nLeft = indices[nCurIdx].first();
            nRight = indices[nCurIdx].second();
            
            // [0,-1] means the entire range.
            
            if(nLeft==0 && nRight == -1)
                nRight = nChildren - 1;
            // Make sure that the indices do not exceed the number of children.
            if(nRight >= nChildren)
                FatalError
                    << "Error when parsing indices: The #" << (nCurIdx+1)
                    << " range ["<< nLeft <<"," << nRight<<"]:\n"
                    << " Cannot update indices bigger than number of children("
                    << nChildren << "). Aborting\n"
                    << exit(FatalError);

            for(label nChildIdx = nLeft; nChildIdx <= nRight; ++nChildIdx)
                children[nChildIdx]->updateNodes(indices,update, nCurIdx+1);
        }
        // Recompute size assuming children are updated.
        if(!children.empty()) {
            size = 0;
            forAll(children, i)
                size += children[i]->size;
        }
    }

    void multiNodeDecomp::nodeMetadata::setLeveledDictionaries(const List<const dictionary*>& dictionaries, label nLevel) {
        // Set the dictionary to this level, and to non-leaf children.
        setDict(*dictionaries[nLevel]);
        forAll(children, i) {
            if(children[i]->nDomains() > 0)
                children[i]->setLeveledDictionaries(dictionaries,nLevel+1);
            
        }
    }

}

// ************************************************************************* //
