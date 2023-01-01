from .model import *


class khipu_diagnosis(Khipu):
    '''Added diagnostic and exploratory functions to khipu class.
    They should be run after khipu.build_khipu().
    '''
    def show_trimming(self):
        '''Show what nodes are removed from input_network. After running self.build_khipu.
        '''
        print("nodes_to_use: ", self.nodes_to_use)
        print("redundant_nodes: ", self.redundant_nodes)

    def build_diagnostic_tree(self, input_network, depth_limit = 10):
        '''Build diagnostic tree using input nodes.

        Updates
        -------
        self.diagnostic_tree as an instance of treelib.Tree.

        Examples (tree here is diagnostic_tree)
        --------
        >>> KP = model.khipu()
        >>> KP.build_tree(big[0], peak_dict, edge_dict)
        >>> 
        >>> KP.tree.show()
        116.0707@138.3
        └── 117.0741@138.5
            ├── 120.084@136.8
            └── 121.0874@136.8

        >>> KP = model.khipu()
        >>> 
        >>> KP.build_tree(big[10], peak_dict, edge_dict)
        []
        >>> 
        >>> KP.tree.show()
        85.0744@110.6
        └── 126.1009@110.6
            └── 127.1043@110.6
                ├── 121.0842@109.9
                ├── 121.0842@110.6
                ├── 122.0876@109.7
                │   ├── 120.0808@109.9
                │   └── 120.0808@110.6
                ├── 122.0876@110.4
                ├── 124.0943@110.8
                └── 125.0976@110.6

        Notes
        -----
        A minimum_spanning_tree will have all necessary patterns to cover full khipu, but not unique pattern.
        '''
        T = nx.minimum_spanning_tree(input_network)
        peak_dict = self.feature_dict
        tree = treelib.Tree()
        root_node = self.root
        tree.create_node(make_peak_tag(peak_dict[root_node]), root_node, data=peak_dict[root_node])

        for node in T[root_node]:
            tree.create_node(make_peak_tag(peak_dict[node]), node, parent=root_node, data=peak_dict[node])

        remaining = [E for E in T.edges() if root_node not in E]
        while remaining and depth_limit > 0:
            tmp = []
            for edge in remaining:
                if edge[0] in tree:
                    node = edge[1]
                    tree.create_node(make_peak_tag(peak_dict[node]), node, parent=edge[0], data=peak_dict[node])
                elif edge[1] in tree:
                    node = edge[0]
                    tree.create_node(make_peak_tag(peak_dict[node]), node, parent=edge[1], data=peak_dict[node])
                else:
                    tmp.append(edge)

            remaining = tmp
            depth_limit -= 1

        if remaining:
            print(remaining)

        return tree

    def build_diagnostic_tree_full(self):
        self.diagnostic_tree = self.build_diagnostic_tree(self.input_network)
        print("Diagnostic tree with all input nodes: ")
        self.diagnostic_tree.show()

    def build_diagnostic_tree_clean(self):
        self.clean_tree = self.build_diagnostic_tree(self.clean_network)
        print("Minimal khipu tree: ")
        self.clean_tree.show()

    def build_khipu_tree(self):
        peak_dict = self.feature_dict
        tree = treelib.Tree()
        root_node = str(round(self.neutral_mass, 4))
        tree.create_node(str(root_node), root_node)
        # tree depth = 2, defined by khipu
        for n in self.khipu_grid.columns:
            pdata = peak_dict.get(n, {})
            plabel = n
            if pdata:
                plabel = make_peak_tag(pdata)
            tree.create_node(plabel, n, parent=root_node, data=pdata)
            for m in self.khipu_grid.index:
                b = self.khipu_grid.loc[m, n]
                if b:
                    tree.create_node(make_peak_tag(peak_dict[b]), b, parent=n, data=peak_dict[b])

        self.khipu_tree = tree
        print("Aligned khipu tree: ")
        self.khipu_tree.show()

    def snap_features_to_grid_by_best_mz(self, expected_grid_mz_values):
        '''Alternative method to create khipu_grid, by best matching each feature to the expected_grid_mz_values.
        This has major problems when there's confusion btw close values, 
        as expected_grid_mz_values can also be shifted by measurement error in root m/z.
        '''
        expected = np.array(expected_grid_mz_values).T
        khipu_grid = pd.DataFrame( np.empty(expected.shape, dtype=np.str),
                            index=self.isotope_index,
                            columns=self.adduct_index,
                            dtype=str)
        for x in self.sorted_mz_peak_ids:
            ii = np.argmin(abs(expected - x[0]))
            khipu_grid.iloc[np.unravel_index(ii, expected.shape)] = x[1]
        return khipu_grid


# -----------------------------------------------------------------------------

def local_read_file(infile, 
                    start_col=3,
                    end_col=6,
                    isotope_search_patterns=isotope_search_patterns, 
                    adduct_search_patterns=adduct_search_patterns,
                    mz_tolerance_ppm=5,
                    rt_tolerance=2, ):
    '''The input feature table must be a tab delimited file, with the first four columns as:
    ID, m/z, retention_time, intensity.
    Example data at '../testdata/full_Feature_table.tsv'.
    '''
    print("Working on file ", infile)
    flist = read_features_from_text(open(infile).read(),
                    id_col=0, mz_col=1, rtime_col=2, intensity_cols=(start_col, end_col), delimiter="\t"
    )
    subnetworks, peak_dict, edge_dict = peaks_to_networks(flist,
                    isotope_search_patterns,
                    adduct_search_patterns,
                    mz_tolerance_ppm,
                    rt_tolerance
    )
    return subnetworks, peak_dict, edge_dict

def khipu_annotate(args):
    '''Automated pre-annotation using Khipu on a feature table.
    args as from parser.parse_args().

    One can follow the steps here to construct a custom workflow,
    using scripting or notebooks.
    '''
    adduct_patterns = adduct_search_patterns
    delta_mz = PROTON
    if args.mode == 'neg':
        adduct_patterns = adduct_search_patterns_neg
        delta_mz = -PROTON

    subnetworks, peak_dict, edge_dict = local_read_file(infile=args.input,
                    start_col=args.start,
                    end_col=args.end,
                    isotope_search_patterns=isotope_search_patterns, 
                    adduct_search_patterns=adduct_patterns,
                    mz_tolerance_ppm=args.ppm,
                    rt_tolerance=args.rtol,
    )
    WV = Weavor(peak_dict, isotope_search_patterns=isotope_search_patterns, 
                adduct_search_patterns=adduct_patterns, 
                mz_tolerance_ppm=args.ppm, mode=args.mode)
    print("\n")
    print("Initial khipu search grid: ")
    print(WV.mzgrid)
    print("\n")

    khipu_list = graphs_to_khipu_list(
        subnetworks, WV, mz_tolerance_ppm=args.ppm,
        )
    khipu_list, all_assigned_peaks = extend_khipu_list(khipu_list, peak_dict, 
                    extended_adducts, mz_tolerance_ppm=args.ppm,
                    rt_tolerance=args.rtol)

    print("\n\n ~~~~~~ Got %d khipus, with %d features ~~~~~~~ \n\n" 
                %(len(khipu_list), len(all_assigned_peaks)))
    empCpds = export_empCpd_khipu_list(khipu_list)

    outfile = 'khipu_test_empricalCompounds.json'
    if args.output:
        outfile = args.output + '.json'
    with open(outfile, 'w', encoding='utf-8') as f:
        json.dump(empCpds, f, ensure_ascii=False, indent=2)

    singleton_ion = list(WV.mzgrid.columns)[0]
    export_khipu_table(outfile.replace(".json", ".tsv"), peak_dict, empCpds, singleton_ion, delta_mz)

def peaklist_to_khipu_list(peaklist, 
                    isotope_search_patterns=isotope_search_patterns, 
                    adduct_search_patterns=adduct_search_patterns,
                    mz_tolerance_ppm=5,
                    rt_tolerance=2, 
                    mode='pos'):
    '''A wrapper function easier for local data analysis.
    peaklist : list of input peaks/features in JSON notion as defined in Khipu or MetDataModel.
    return khipu_list, all_assigned_peaks
    '''
    subnetworks, peak_dict, edge_dict = peaks_to_networks(peaklist,
                    isotope_search_patterns,
                    adduct_search_patterns,
                    mz_tolerance_ppm,
                    rt_tolerance
    )
    WV = Weavor(peak_dict, isotope_search_patterns=isotope_search_patterns, 
                adduct_search_patterns=adduct_search_patterns, 
                mz_tolerance_ppm=mz_tolerance_ppm, 
                mode=mode)
    print("\n")
    print("Initial khipu search grid: ")
    print(WV.mzgrid)
    print("\n")
    khipu_list = graphs_to_khipu_list(
        subnetworks, WV, mz_tolerance_ppm=mz_tolerance_ppm,
        )
    khipu_list, all_assigned_peaks = extend_khipu_list(khipu_list, peak_dict, 
                    extended_adducts, mz_tolerance_ppm=mz_tolerance_ppm,
                    rt_tolerance=rt_tolerance)
    return khipu_list, all_assigned_peaks

def graphs_to_khipu_list(subnetworks, WeavorInstance, mz_tolerance_ppm):
    '''Generate full khipu_list from subnetworks, 
    including iterative khipus based on features pruned out of initial subnetwork.
    '''
    khipu_list = []
    for g in subnetworks:
        if g.number_of_edges() > 0:
            KP = Khipu(g)
            KP.build_khipu(WeavorInstance, mz_tolerance_ppm)
            khipu_list.append(KP)
            while KP.redundant_nodes and KP.pruned_network and KP.pruned_network.edges():
                # make sure pruned_network is connected
                more_subnets = [KP.pruned_network.subgraph(c).copy() 
                                        for c in nx.connected_components(KP.pruned_network)]
                for _G in more_subnets:
                    KP = Khipu(_G, isotope_search_patterns, adduct_search_patterns)
                    KP.build_khipu(WeavorInstance, mz_tolerance_ppm)
                    khipu_list.append(KP)

    # assign IDs
    ii = 0
    for KP in khipu_list:
        base_mz = str(round(KP.neutral_mass, 4))
        ii += 1
        KP.id = 'kp' + str(ii) + "_" + str(base_mz)

    return khipu_list

def extend_khipu_list(khipu_list, peak_dict, adduct_search_patterns_extended, mz_tolerance_ppm=5, rt_tolerance=2):
    '''Update khipus by extended adduct search.
    Returns updated khipu_list and list_assigned_peaks.
    '''
    list_assigned_peaks = []
    for KP in khipu_list:
        list_assigned_peaks += KP.nodes_to_use
    unassigned_peaks = [v for x,v in peak_dict.items() if x not in list_assigned_peaks]
    mztree = build_centurion_tree(unassigned_peaks)
    for KP in khipu_list:
        added_peaks = KP.extended_search(mztree, 
                            adduct_search_patterns_extended,  mz_tolerance_ppm, rt_tolerance)
        list_assigned_peaks += added_peaks

    return khipu_list, set(list_assigned_peaks)

def export_json_khipu_list(khipu_list):
    J = []
    for KP in khipu_list:
        J.append(KP.export_json())
    return J

def export_empCpd_khipu_list(khipu_list):
    '''Export all khipus in khipu_list to a list of empirical compounds, which is JSON compatible.
    Wrapper of khipu.format_to_epds().
    A small number of features are "undetermined", as they came due to initial edges but violate DAG rules. 
    They are sent off for a new khipu.
    '''
    J = []
    for KP in khipu_list:
        J.append(KP.format_to_epds(id=KP.id))
    return J

def export_khipu_table(outfile, peak_dict, json_khipu_list, singleton_ion='M+H+', delta_mz=1.0073):
    '''Write tab delimited file of khipu_list.
    delta_mz is +proton if pos mode, -protone if neg mode, for singletons.
    '''
    header = ['id', 'mz', 'rtime', 'interim_id', 'neutral_formula_mass', 'isotope', 'modification', 'ion_relation']
    s = '\t'.join(header).replace('interim_id', 'empCpd') + '\n'
    in_khipu_list = []
    for J in json_khipu_list:
        for F in J["MS1_pseudo_Spectra"]:
            in_khipu_list.append(F['id'])
            s += '\t'.join([
                F['id'], str(F['mz']), str(F['rtime']), J['interim_id'], str(J['neutral_formula_mass']), 
                F.get('isotope', ''), F.get('modification', ''), F['ion_relation']
            ]) + '\n'
    # other peaks
    for k,F in peak_dict.items(): 
        if k not in in_khipu_list:
            s += '\t'.join([
                F['id'], str(F['mz']), str(F['rtime']), 'singleton', str(round(F['mz']-delta_mz,4)), 
                'M0', singleton_ion, ''
            ]) + '\n'

    with open(outfile, 'w') as f:
        f.write(s)
