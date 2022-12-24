
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


class not_use:
    def build_grid_abstracted(self, abstracted_adduct_edges, root_branch):
        '''Compute grid structure and expected_grid_mz_values. 
        This uses DAG traversal but more complex than necessary, and does not find neutral mass.

        Parameters
        ----------
        abstracted_adduct_edges : directed edges that connect B-nodes based on adduct m/z difference not isotopes.
            B-nodes are abstracted branches from isotopes during grid building.

        Returns
        -------
        indexed_adducts : a list containing ordered adducts (some are abstracted branches)
        adduct_index_labels : annotation companion to adduct_index, list of same node order
        expected_grid_mz_values : calculated m/z values for khipu grid based on self.root, indexed_adducts 
            and self.isotope_search_patterns.

        Notes
        -----
        To compute m/z from abstracted_adduct_edges, we need to determine the nodes for khipu trunk and their order.
        Build a DiGraph which should be DAG. Use nx.dfs_edges to walk the DAG.
        Grid is computed on the trunk (from walking DAG) and isotope_search_patterns.

        '''
        root_mz = self.feature_dict[self.root]['mz']
        adduct_mz_dict, mz_modification_dict, node_mz_dict = {}, {}, {root_branch: root_mz}
        dict_node_label = {}
        for a in self.adduct_search_patterns:
            adduct_mz_dict[a[1]] = a[0]
        for e in abstracted_adduct_edges:
            mz_modification_dict[make_edge_tag((e[0], e[1]))] = adduct_mz_dict[e[2]['tag']]
            # e[1] + 
            dict_node_label[e[1]] = e[0] + '+' + e[2]['tag']

        DG = nx.DiGraph(abstracted_adduct_edges)
        # walk the graph through all nodes, and build trunk in that order
        indexed_adducts, trunk_mzlist = [root_branch, ], [root_mz, ]

        try:
            for e in nx.dfs_edges(DG, source=root_branch):
                indexed_adducts.append(e[1])
                target_mz = node_mz_dict[e[0]] + mz_modification_dict[make_edge_tag((e[0], e[1]))]
                # e[0] is guanranteed to be in node_mz_dict due to the tree walking order
                node_mz_dict[e[1]] = target_mz
                trunk_mzlist.append(target_mz)

        except KeyError:
            print("DAG violation: ", root_mz, self.root)
            self.exceptions.append("DAG violation")
            # print(DG.edges())
            G = nx.Graph(abstracted_adduct_edges)
            indexed_adducts, trunk_mzlist = [root_branch, ], [root_mz, ]
            for e in nx.bfs_edges(G, source=root_branch):
                indexed_adducts.append(e[1])
                target_mz = node_mz_dict[e[0]] + mz_modification_dict[make_edge_tag((e[0], e[1]))]
                node_mz_dict[e[1]] = target_mz
                trunk_mzlist.append(target_mz)

        adduct_index_labels = [dict_node_label.get(x, x) for x in indexed_adducts]
        expected_grid_mz_values = []
        for A in trunk_mzlist:
            expected_grid_mz_values.append(
                [A] + [A+x[0] for x in self.isotope_search_patterns]
            ) 

        return indexed_adducts, adduct_index_labels, expected_grid_mz_values


    def snap_features_to_grid(self, branch_dict, expected_grid_mz_values):
        '''Create khipu_grid. To snap each feature to the expected_grid_mz_values.

        Parameters
        ----------
        branch_dict : dictionary, branch ID to member features/nodes.
        expected_grid_mz_values : m/z values in a list of branches. 
            List of lists, _N x _M, in same order as self.adduct_index.

        Updates
        -------
        self.annotation_dict : the isotope and adduct notions for each feature.
        self.redundant_nodes : features that do not fit DAG are sent to redundant_nodes.

        Returns
        -------
        khipu_grid : pd.DataFrame, adducts as cols (trunk) and isotopes as rows (branches).

        Notes
        -----
        DataFrame is better for khipu_grid, because easier to cast data types and keep matrix format.
        Nested list or numpy array is harder to ensure correct matrix format.
            The algorithm here has to respect established edges to avoid confusion. 
        self.adduct_index and branch_dict define memberships in each adduct group.
        The method below is very error prone:
            for x in self.sorted_mz_peak_ids:
                ii = np.argmin(abs(expected - x[0]))
                khipu_grid.iloc[np.unravel_index(ii, expected.shape)] = x[1]

        If isotope_search_patterns and adduct_search_patterns are not used propersly,
        each pattern may occur more than once. Not all edges in input network belong to this Khipu.
        Tree search methods don't always cover all edges, not a good choice.
        If minimum_spanning_tree is calculated from self.input_network, the isotopic branch may not be fully connected.

        '''
        _M, _N = len(self.isotope_index), len(self.adduct_index)
        khipu_grid = pd.DataFrame( np.empty((_M, _N), dtype=np.str),
                            index=self.isotope_index,
                            columns=self.adduct_index_labels ,                  # adduct_index
                            dtype=str)
        for jj in range(_N):
            # get cooridnate for each matched isotope
            if self.adduct_index[jj] in branch_dict:
                isotopes = branch_dict[self.adduct_index[jj]]
            else:
                isotopes = [self.adduct_index[jj]]                      # single feature
            for F in isotopes:
                ii = np.argmin([abs(x-self.feature_dict[F]['mz']) for x in expected_grid_mz_values[jj]])
                khipu_grid.iloc[ii, jj] = F
                self.annotation_dict[F] = (self.isotope_index[ii], self.adduct_index_labels[jj])
        
        for n in self.nodes_to_use:
            if n not in self.annotation_dict:
                self.redundant_nodes.append(n)

        return khipu_grid

