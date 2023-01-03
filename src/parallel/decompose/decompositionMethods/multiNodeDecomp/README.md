# New Multi-Level Decomposition
The multi-node decomposition is an extension of the existing multi-level decomposition. It supports the syntax of the current multi-level decomposition, but allows to change the decomposition tree as you wish. For example, you may split into unbalanced nodes, set the weights of some nodes to be bigger than others, or perhaps use a different decomposition method for some nodes.
You may set up the decomposition in two ways:
1. Using a domains list and a default method:
    ```
    numberOfSubdomains 8;
    multiNodeCoeffs {
        domains (2 4);
        method metis;
    }
    ```
2. Using a dictionary for each level:
    ```
    numberOfSubdomains 8;
    multiLevelCoeffs {
        nodes {
            numberOfSubdomains 2;
            method metis;
        }
        cores {
            numberOfSubdomains 4;
            method scotch;
        }
    }
    ```
    Note that if the total number of subdomains does not match the product of the number of subdomains at each level, but a default method is provided, a new level will be inferred in order to match the total number of subdomains.

This creates a "decomposition tree" - for example, the dictionaries above create a tree, where the root has two children, and each child has four children (who are the leaves of the tree). Every leaf in the tree is a subdomain in the final decomposition.
After setting up the decomposition, we may edit specific nodes or ranges of nodes. For example, suppose we want to split into two nodes, the first one having four subdomains and the second having eight subdomains. We can use the above dictionaries, and then use:
```
domains[1] (8);
```
The squared brackets indicate which nodes in the tree should we edit - We want the second child of the root (the indexing starts from zero). If we wanted to change the first two children of the third child of the root, we would write:
```
domains[2][0-1] (8);
```

Note that the total number of subdomains must match the number of subdomains declared after all modifications. In addition, note that the decomposition into two nodes will be done as if they were of the same size, hence the first four subdomains will be bigger than the other eight. In order to fix this, we may:
1. Change the weight of the second node into twice the weight:
```
weight[1] 2;
```
2. Set the weights initialization into relative - this will cause the weights of the children to first be computed by the amount of leaves in their subtree. Note that this updates the whole subtree initialization, but using the `weight` parameter, we can override this initialization.
```
weightsInitialization[1] relative;
```


We may also set a special method dictionary that decomposes differently for some nodes:
```
method[2-4] {
    numberOfSubdomains 4;
    method metis;
    coeffs {
        ...
    }
}
```
