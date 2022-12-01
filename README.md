# khipu as a tree data structure to model a chemical compound in mass spectrometry data

This applies to regular LC-MS data, but also enables easy analysis of isotope tracing and chemical labeling data.


## Algorithmic consideration: the depth of trees: 
Relation btw different number of labels? M3 is not parent of M4 - they are children of root and at the same level.
How many levels of adducts? This should be 1 because adduct formation is one-time event in ESI.
Thus, annoTree should have 2-levels: isotopes and adducts. 
As we separate search patterns from search functions, 
one can write separate logic to design search patterns and further process results.

For chemical derivatization, the data can be organized as tree pairs, without changing the trees internally.
This is because derivatization is a reaction that occurs before LC-MS, and
LC-MS measures whatever compounds that are present in samples.

Demonstrating data and concepts:

    >>> node2tree = construct_isotopic_trees(plist)
    Found 873 isotopic pairs, 703 trees and 170 in branches.
    >>> node2tree['F795'].root
    'F649'
    >>> node2tree['F795'].nodes
    {'F649': Node(tag=179.1123@172.9, identifier=F649, data=None), 'F755': Node(tag=180.1156@172.9, identifier=F755, data=None), 
    'F794': Node(tag=181.1198@171.7, identifier=F794, data=None), 'F795': Node(tag=181.1198@173.3, identifier=F795, data=None), 
    'F839': Node(tag=182.1233@171.7, identifier=F839, data=None)}
    >>> node2tree['F795'].show()
    179.1123@172.9   #  (M0)
    └── 180.1156@172.9
        ├── 181.1198@171.7
        │   └── 182.1233@171.7    # (should be level 1 not lower)
        └── 181.1198@173.3        # problematic assignment

## Algorithmic consideration: tree search vs network partitioning.
We can 
a) search all isotopic trees and merge the adducts in two steps, or 
b) search all isotopic pairs and all adduct pairs, form the network and partition the network.
The a) option may have problem if the merging only considers M0, 
which may not be the most abundant ion in labeling experiments. 
I.e. M0 may not be measured for a compound and the adduct relation of M0s is broken.
If we have to search beyond M0 for adduct relations, a) offers no advantage.
Therefore, b) is the more generic solution.

## Algorithmic consideration: assigning nodes to khipu grid.
We can 

    >>> big[0].edges()
    [('F1478', 'F171'), ('F1478', 'F114'), ('F1478', 'F10'), ('F171', 'F114'), ('F171', 'F10'), ('F114', 'F10')]
    >>> T = nx.minimum_spanning_tree(big[0])
    >>> 
    >>> T.edges()
    [('F1478', 'F171'), ('F1478', 'F114'), ('F1478', 'F10')]
    >>> 


## Implementation
1. Get all pairs of features that match isotope patterns and adduct patterns. The patterns are provided as options and users can use their own patterns. The default match in LC-MS is based on m/z and retention time windows, while preprocessing software (asari) can require overlap scans.

2. Connect all pattern matched feature pairs to an overall network, which is further partitioned into connected subnetworks.

3. Each subnetwork becomes a khipu instance, after rule inspection and optimization to resolve ambiguity. A khipu is essentially an empirical compound that is used for downstream annotation and analysis.

4. Visualization within khipu library is limited to tree and matrix printing in plain text. Plot functions of isotopologues and adducts will be in plotly and notebooks, not in this library.


## The khipus are used as empirical compounds in our computational libraries (https://github.com/shuzhao-li/metDataModel), which can be annotated by customized methods, e.g.
- in-house compound library
- targeted pathway search
- public databases
- linking tandem mass spec data


