
#
# -----------------------------------------------------------------------------
#
# Not currently used
#

def construct_deep_isotopic_trees(peak_list, 
                    search_patterns = [ (1.003355, '13C/12C', (0, 0.8)),
                                        (2.00671, '13C/12C*2', (0, 0.8)),
                                        (3.010065, '13C/12C*3', (0, 0.8)),
                                        (4.01342, '13C/12C*4', (0, 0.8)),
                                        (5.016775, '13C/12C*5', (0, 0.8)),
                                        (6.02013, '13C/12C*6', (0, 0.8)),], 
                    mz_tolerance_ppm=5, 
                    isotope_rt_tolerance=2, 
                    check_isotope_ratio = False,
                    tree_depth_limit=10):
    '''
    Make deep isotopic trees from peak_list, where M3 is treated as child of M2.
    This is not preferred. 
    tree_depth_limit: limit of tree depth, i.e. steps on a branch.
    For other parameters, see get_isotopic_pairs.
    Return node2tree, a dict of peak to tree mapping.
    No peak is assigned to more than one trees. But some close peaks can't be distinguished here.
    search_patterns: Need to consider all possible number of 13C labels. 
                     If *6 exists, it's not granted that *4 exists. So we can't rely on pairwise to connect all pairs.
                     The third item, (0, 0.8) here, is an option to constrain ratios but not used in this function.
    Example
    =======
    >>> node2tree['F515'].show()
    173.0922@172.9
    └── 174.0956@171.5
        ├── 175.0991@172.4
        │   ├── 176.1031@171.7
        │   │   ├── 177.1065@171.7
        │   │   └── 177.1065@174.1
        │   └── 176.1031@174.1
        └── 175.0991@173.1
    >>> node2tree['F796'].show()
    180.1244@171.2
    └── 181.1279@170.1
        └── 182.1313@170.1
    '''
    peak_dict = make_peak_dict(peak_list)
    mztree = build_centurion_tree(peak_list)
    isotopologues = get_isotopic_pairs(peak_list, mztree, search_patterns, mz_tolerance_ppm, 
                    isotope_rt_tolerance, check_isotope_ratio)
    # isotopologues format: [ ((195, 'anchor'), (206, '13C/12C')), ...]
    # build trees
    annoTrees, branches = [], []
    all_target_nodes = set([x[1][0] for x in isotopologues])
    for pair in isotopologues:
        if pair[0][0] not in all_target_nodes:   # if source_node appears in any of the target_nodes, it's not a root
            tree = treelib.Tree()
            tree.create_node(make_peak_tag(peak_dict[pair[0][0]]), pair[0][0], data=pair[0][1])
            tree.create_node(make_peak_tag(peak_dict[pair[1][0]]), pair[1][0], parent=pair[0][0], data=pair[1][1])
            annoTrees.append(tree)
        else:
            branches.append(pair)
    
    print("Found %d isotopic pairs, %d trees and %d in branches." %(len(isotopologues), len(annoTrees), len(branches)))
    
    node2tree = {}
    for tree in annoTrees:
        for node in tree.nodes:
            node2tree[node] = tree

    # do branches now
    remaining = []
    for pair in branches:          
        if pair[0][0] in node2tree:                  # pair[0] already in a tree
            try:
                this_tree = node2tree[pair[0][0]]
                this_tree.create_node( make_peak_tag(peak_dict[pair[1][0]]), pair[1][0], parent=pair[0][0], data=pair[1][1] )
                node2tree[pair[1][0]] = this_tree
            except treelib.exceptions.DuplicatedNodeIdError:
                # print("already included ", pair)
                pass                 # pair already in a tree
        else:
            remaining.append(pair)

    steps = 0
    while remaining and steps < tree_depth_limit:
        tmp = []
        for pair in remaining:          
            if pair[0][0] in node2tree: 
                try:
                    this_tree = node2tree[pair[0][0]]
                    this_tree.create_node( make_peak_tag(peak_dict[pair[1][0]]), pair[1][0], parent=pair[0][0], data=pair[1][1] )
                    node2tree[pair[1][0]] = this_tree
                except treelib.exceptions.DuplicatedNodeIdError:
                    # print("already included ", pair)
                    pass               # pair already in a tree
            else:
                tmp.append(pair)
        remaining = tmp
        steps += 1

    return node2tree

def merge_trees_by_modifications(trees, list_peaks, search_patterns, rt_verify_function,
                    mz_tolerance_ppm=5, rt_tolerance=10):
    '''
    Merge trees that are from same compound but with  modifications. 
    rt_verify_function applies rules on retention time.
    '''
    peak_dict = make_peak_dict(list_peaks)
    tree_dict = {}
    for tree in trees: tree_dict[tree.root] = tree
    print("Merging adducts on %d trees..." %len(trees))
    matched_pairs = []
    relevant_peak_centTree = build_centurion_tree( [peak_dict[p] for p in tree_dict.keys()] )
    for p in tree_dict:
        P1 = peak_dict[p]
        matched = [  ] 
        for _pair in search_patterns:
            (mass_difference, relation) = _pair[:2]
            tmp = find_all_matches_centurion_indexed_list(P1['mz'] + mass_difference, 
                                                    relevant_peak_centTree, mz_tolerance_ppm)
            for P2 in tmp:
                if rt_verify_function(P1, P2, rt_tolerance):
                    matched.append( (P1['id'], P2['id'], relation) )
        matched_pairs += matched

    list_primary = set([x[0] for x in matched_pairs])
    list_adducts = set([x[1] for x in matched_pairs])
    overlap = list_primary.intersection(list_adducts)
    union = list_primary.union(list_adducts)
    if overlap: print("Unresolved multiple relationships: ", overlap)

    good_trees = [tree_dict[x] for x in tree_dict if x not in union]
    for P in matched_pairs:
        p1, p2 = P[:2]
        tree_dict[p2].get_node(p2).data = P[2]      # add relation note
        try:
            tree_dict[p1].paste(p1, tree_dict[p2])
        except ValueError:                  # likely Duplicated nodes
            # print(tree_dict[p1].show(), tree_dict[p2].show())
            for node in tree_dict[p2].nodes:
                if node not in tree_dict[p1].nodes:
                    tree_dict[p1].add_node(tree_dict[p2].nodes[node], parent=tree_dict[p1].root)
        good_trees.append(tree_dict[p1])

    print("Got %d merged trees." %len(good_trees))
    return good_trees


def merge_trees_by_insrc_modifications(trees, list_peaks,
                    search_patterns = [(1.0078, 'H'), (21.9820, 'Na/H'), (41.026549, 'Acetonitrile')], 
                    mz_tolerance_ppm=5, rt_tolerance=10):
    '''
    Merge trees that are from same compound but with diff adducts or neutral loss,
    with user supplied search_patterns (can be from search.common_adducts).
    These in-source modifications must have same retention time.
    '''
    trees = merge_trees_by_modifications(trees, list_peaks,
                    search_patterns, rt_verify_function=rt_matched_by_tolerance,
                    mz_tolerance_ppm=mz_tolerance_ppm, rt_tolerance=rt_tolerance)
    return add_data_to_tag(trees)

def merge_trees_by_derivatization(trees, list_peaks,
                    search_patterns = [(161.08407, 'DmPA'), (233.05105, 'Dens'), (247.07793, 'DnsHz')], 
                    mz_tolerance_ppm=5, ):
    '''
    Merge trees that are from same compound but with chemical labeling/derivatization,
    which usually leads to greater retention time.
    '''
    trees = merge_trees_by_modifications(trees, list_peaks,
                    search_patterns, rt_verify_function=rt_compared_by_values,
                    mz_tolerance_ppm=mz_tolerance_ppm, rt_tolerance=None)
    # return add_data_to_tag(trees)
    return trees

