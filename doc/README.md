# Additional documentation for developers

## Why tree?
Grouping ions in mass spectrometry data to their origin compound is a key step of data annotation.
It's a common practice to search the data for signatures of mass differences, 
then combine overlap signatures as a network of ions (also called features or peaks).
Ideally, one wants to infer the neutral formula of the compound, 
and annotate how each ion is related to the origin compound. 
In such case, each ion should have only one unique "parent" in the mass difference network, 
and the network should have one origin compound as root. 
This becomes a tree structure.

In theory, multiple ions are measured in mass spectrometry for one origin compound,
because of two factors:
1) isotopes exists for atoms, which have difference mass values; and
2) modifications occur during the mass spectrometry analysis, including formation of adducts, ionization, in-source fragmentation and neutral loss.

Each ion is a result of combination of these two factors, and can be attributed to the origin compound.
Therefore, the ions can be adequately represented as a "leaf" on a tree, root being the origin compound.

Once we move beyond generic networks and use trees to annotate ions,
it is much easier to optimize algorithms and exchange annotations.

## Tree and depth
After we get paired relationship between ions by mass differences, a tree may look like this:

    179.1123@172.9   #  (M0)
    └── 180.1156@172.9
        ├── 181.1198@171.7
        │   └── 182.1233@171.7    # (should be level 1 not lower)
        └── 181.1198@173.3        # problematic assignment

Here, each ion is denoted as m/z@retention_time. 
The above tree has a problem of tree depth.
For isotopes, M3 is not parent of M4 - they are children of root and at the same level.
For adducts, they should be at the same level because adduct formation is one-time event in ESI.
Thus, the annotation tree should have 2-levels: isotopes and adducts (and other modifications). 

The above tree has an additional problem, as we have two ions at 181.1198, 
one of which could be from a different origin compound.

## Algorithmic consideration: tree search vs network partitioning
Let's consider two senarios of constructing a khipu per origin compound:
a) search all isotopic trees and merge the adducts in an additional step, or 
b) search all isotopic pairs and all adduct pairs, form the network and partition the network.

The a) option may have problem if the merging only considers M0, 
which may not be the most abundant ion in isotope labeling experiments. 
I.e. M0 may not be measured for a compound and the adduct relation of M0s is broken.
If we have to search beyond M0 for adduct relations, a) offers no advantage.
Therefore, b) is the more generic solution.

The network formed by m/z pairs is redundant, especially when multiple isotopes exist. 
An example is that only half of edges are needed to connect the ions:

    >>> big[0].edges()
    [('F1478', 'F171'), ('F1478', 'F114'), ('F1478', 'F10'), ('F171', 'F114'), ('F171', 'F10'), ('F114', 'F10')]
    >>> T = nx.minimum_spanning_tree(big[0])
    >>> 
    >>> T.edges()
    [('F1478', 'F171'), ('F1478', 'F114'), ('F1478', 'F10')]
 
In theory, one can find the ion of lowest m/z, calculate all the combinations of isotopes and adducts, 
and match all ions to the calculated combinations (a grid). 
However, it will be problematic because each measurement has an error, 
and this method imposes the error of one ion to others, sufficient to cause wrong matches 
when m/z patterns are close. E.g. the m/z difference of 13C/12C is 1.003355 and that of a proton is 1.0073.
For high m/z values, their difference is sensitive to measurement errors. 


## Khipus are used as empirical compounds 
- Empirical compound is defined at https://github.com/shuzhao-li/metDataModel
- Khipu defines the MS1 pseudo spectra for empirical compounds
- Empirical compounds are designed to include additional annotations, including from tandem mass spec data and database searches
- Khipu is used in JMS and asari (Trackable and scalable Python program for high-resolution LC-MS metabolomics data preprocessing https://github.com/shuzhao-li/asari). 
