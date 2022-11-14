'''
Comprehensive construction of empCpds via generic tree structures.
Each annoTree = empCpd["MS1_pseudo_Spectra"].

1. construct all isotopologues into a list of trees
2. attach in-src modifications (adducts) to nodes in each tree
3. for peaks without isotopologues, add to new trees and attach adducts

This applies to regular LC-MS data, but also enables easy analysis of isotope tracing and chemical labeling data.

On the depth of trees: 
Relation btw different number of labels? M3 is not parent of M4 - they are children of root and at the same level.
How many levels of adducts? This should be 1 because adduct formation is one-time event in ESI.
Thus, annoTree should have 2-levels: isotopes and adducts. 
As we separate search patterns from search functions, 
one can write separate logic to design search patterns and further process results.

For chemical derivatization, the data can be organized as tree pairs, without changing the trees internally.
This is because derivatization is a reaction that occurs before LC-MS, and
LC-MS measures whatever compounds that are present in samples.


The trees (empCpds) can be annotated by customized methods, e.g.
- in-house compound library
- targeted pathway search
- public databases
- linking tandem mass spec data

Algorithmic consideration: tree search vs network partitioning.
We can 
a) search all isotopic trees and merge the adducts in two steps, or 
b) search all isotopic pairs and all adduct pairs, form the network and partition the network.
The a) option may have problem if the merging only considers M0, 
which may not be the most abundant ion in labeling experiments. 
I.e. M0 may not be measured for a compound and the adduct relation of M0s is broken.
If we have to search beyond M0 for adduct relations, a) offers no advantage.
Therefore, b) is the more generic solution.



# to-dos: 
# To add plot functions of isotopologues



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



Trees use peak IDs, but full peak data are in json peak lists.


# Above only includes isotopic peaks; will add peaks without C13 counterparts into empCpd trees

# from itertools import combinations
'''


import json
import treelib




def network_to_trees(subnetwork, 
                isotope_search_patterns, adduct_search_patterns, peak_dict, edge_dict):
    '''
    Convert a network of two types of edges, isotopic and adduct, to empTree instances of unique nodes.
    Use the grid to enforce unique ion in each position.
    
    to 
    
    '''
    trees = []




    for pair in isotopologues:
        if pair[0][0] not in all_target_nodes:   # if source_node appears in any of the target_nodes, it's not a root
            tree = treelib.Tree()
            tree.create_node(make_peak_tag(peak_dict[pair[0][0]]), pair[0][0], data=pair[0][1])
            tree.create_node(make_peak_tag(peak_dict[pair[1][0]]), pair[1][0], parent=pair[0][0], data=pair[1][1])
            annoTrees.append(tree)
        else:
            branches.append(pair)
    


    return trees













class empTree:
    '''
    Tree representation of an empirical compound.
    The basic functions are in mass2chem.annotree, which uses treelib for tree computing.

    An empTree follows the matrix notion of a family of peaks: a tree of 2 levels - isotopes and adducts.

    It's not always easy or feasible to compute unambiguous tree for an empirical compound, e.g. one may need to consolidate 
        118.0652@109.9
        ├── 119.0685@109.2
        └── 159.0917@109.9
            └── 160.0952@110.4
        159.0917@110.4
        └── 160.0952@110.4

    This class can be used future research to improve the tree constructions. 


    Features of low abundance do not have detectable isotopes, but can have multiple adducts.
    This function is used after isotopic trees are searched exhaustively. Tree depth here is limited to 1.
    Return node2tree, a dict of peak to tree mapping.


    '''
    def __init__(self):
        '''
        
        '''
        self.tree = None
        self.uniqe_nodes = []
        self.edges = []
        self.is_13C_verified = False
        self.is_singleton = False



    def organize(tree):
        '''
        This is the arborist of the tree
        To clean up a tree, by resovling "sick branches", i.e. misfit redundant features.
        Each m/z should have only one feature to include in a tree.
        Fix tags also.
        Return a list of trees, since a cut branch could become a new tree.

        '''
        pass


    def format_to_epds(self):
        '''
        Add to 'ion_relation' - 'isotope': [M0, C13*M3..], 'modification': ['-CO2', ...]






        '''
        pass





