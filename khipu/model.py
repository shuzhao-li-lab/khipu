'''
Comprehensive construction of empCpds via generic tree structures (khipu).
Each khipu = empCpd["MS1_pseudo_Spectra"].

To-Dos:

- cli function, pos/neg
- automate ext search

'''

import numpy as np
import pandas as pd

try:
    import matplotlib.pyplot as plt
except:
    print("  matplotlib ImportError.")

from .utils import *

class khipu:
    '''2-tier tree representation of an empirical compound.
    An empTree follows the matrix notion of a family of peaks: a tree of 2 levels - isotopes and adducts.
    A true root of khipu is the compound in neutral formula.
    A peudo root of khipu is the ion of lowest m/z.
    The representative adduct should have the most abundance.
    The matrix notion of a khipu is a DataFrame: columns are Adduct labels; rows isotopes.
    '''
    def __init__(self, subnetwork, isotope_search_patterns, adduct_search_patterns):
        '''Do not include negative m/z difference in this class, as the node of lowest m/z is used as root.
        Neutral loss and fragments can be searched after khipu is established.
        Input subnetwork is undirected graph, but khipu will require directed graph/tree, by increasing m/z.
        subnetwork : undirected graph. Example edges:
            ('F3533', 'F20', {'type': 'modification', 'tag': 'Na/H'}), 
            ('F195', 'F20', {'type': 'modification', 'tag': 'Acetonitrile'}), 
            ('F20', 'F53', {'type': 'isotope', 'tag': '13C/12C'}), 
        
        self.is_13C_verified = False
        self.is_singleton = False
        '''
        self.input_network = subnetwork
        self.isotope_search_patterns = isotope_search_patterns
        self.adduct_search_patterns = adduct_search_patterns

        self.nodes_to_use = []
        self.annotation_dict = {}
        self.redundant_nodes = []
        self.pruned_network = None

        self.adduct_index = []
        self.isotope_index = ['M0'] + [x[1] for x in isotope_search_patterns]
        self.khipu_grid = []                # DataFrame of N x M, M number of isotopes, N number of adducts
                                            # M is explicitly specified during search
                                            # N is dynamically determined on each khipu


    def build_khipu(self, peak_dict, mz_tolerance_ppm=5, check_fitness=False):
        '''Convert a network of two types of edges, isotopic and adduct, to khipu instances of unique nodes.
        Use the grid to enforce unique ion in each position, as the initial network can contain 
        erroneously redundant leaves/nodes.

        check_fitness : if True, chech if replacing with any redundant nodes improves fitness. 
            This is a placeholder as a good fitness function is not easy to implement and it slows down computing.

        Updates
        -------
        self.khipu_grid : pd.DataFrame
        self.pruned_network : extra features and edges that are not fit in this khipu
        '''
        self.feature_dict, self.mzstr_dict = self.get_feature_dict(peak_dict, mz_tolerance_ppm)
        if self.input_network.number_of_edges() == 1:
            self.build_simple_pair_grid(peak_dict)
            self.clean_network = self.input_network
        else:
            self.clean()
            self.clean_network = self.input_network.subgraph(self.nodes_to_use)
            isotopic_edges, adduct_edges = [], []
            for e in self.clean_network.edges(data=True):
                if e[2]['type'] == 'isotope':
                    isotopic_edges.append(e)
                else:                   # no other type of connections are allowed 
                    adduct_edges.append(e)
            
            if isotopic_edges and adduct_edges:
                self.build_full_grid(isotopic_edges, adduct_edges)
            elif isotopic_edges :       # branch only A1 adduct
                self.build_branch_only_grid()
            else:                       # trunk only, no isotopes
                self.build_trunk_only_grid(adduct_edges)
            
            if check_fitness:
                self.select_by_fitness()

        if self.redundant_nodes:
            self.pruned_network = self.input_network.subgraph(self.redundant_nodes)

    def build_simple_pair_grid(self, peak_dict):
        '''A khipu of only two nodes (one edge) does not need to go through full grid process,
        and is formatted here. The full grid process would work for these simple cases, but less efficient.
        '''
        e = self.input_network.edges(data=True)[0]
        if peak_dict[e[0]]['mz'] > peak_dict[e[1]]['mz']:
            e = (e[1], e[0], e[2])
        self.root = e[0]
        if e[2]['type'] == 'isotope':
            self.khipu_grid = pd.DataFrame( {'A1': [e[0], e[1]]},
                            index=['M0', e[2]['tag']], 
                            dtype=str)
            self.annotation_dict[e[0]] = ('M0', 'A1')
            self.annotation_dict[e[1]] = (e[2]['tag'], 'A1')
        else:           # adduct
            self.khipu_grid = pd.DataFrame( {'A1': [e[0]], e[2]['tag']: [e[1]]},
                            index=['M0',], 
                            dtype=str)
            self.annotation_dict[e[0]] = ('M0', 'A1')
            self.annotation_dict[e[1]] = ('M0', e[2]['tag'])

    def clean(self):
        '''Clean up the input subnetwork, only using unique features to build a khipu frame.
        Redundant features are kept aside. The leftover features will be sent off to build new khipus.

        Notes
        -----
        Note: mzstr_dict can be problematic in data that were not processed well, 
        because minor potential shift can confuse ion relations.
        This clean() may cut some initial edges.
        '''
        self.median_rtime = np.median([self.feature_dict[n]['rtime'] for n in self.input_network])
        for k,v in self.mzstr_dict.items():
            # v as list of node IDs
            if len(v) == 1:
                self.nodes_to_use.append(v[0])
            else:
                sorted_by_rtime_match = sorted(
                    [(abs(self.feature_dict[n]['rtime'] - self.median_rtime), n) for n in v]
                )
                self.nodes_to_use.append(sorted_by_rtime_match[0][1])
                self.redundant_nodes += [n[1] for n in sorted_by_rtime_match[1:]]

        self.sorted_mz_peak_ids = self.sort_nodes_by_mz()
        self.root = self.sorted_mz_peak_ids[0][1]

    def get_feature_dict(self, peak_dict, mz_tolerance_ppm):
        '''Index all input features; establish str identifier for features of same/close m/z values.
        Base on asari mass track IDs; keep unique m/z only.
        It's more efficient to use, since feature_dict is much smaller than peak_dict.

        Parameters
        ----------
        peak_dict : dict of peaks/features indexed by IDs. Must have fields 
                    'id', 'mz', 'rtime', 'representative_intensity'.
        mz_tolerance_ppm : ppm tolerance in examining m/z groups.

        Returns
        -------
        feature_dict : feature_dict with 'parent_masstrack_id'.
        mzstr_dict : dict indexed by parent_masstrack_id, e.g. {'trackx': [f1, f2], ...}
        '''
        feature_dict, mzstr_dict = {}, {}
        for n in self.input_network.nodes():
            feature_dict[n] = peak_dict[n]
        if 'parent_masstrack_id' not in feature_dict[n]:
            feature_dict = assign_masstrack_ids_in_khipu(feature_dict, mz_tolerance_ppm)

        for n,v in feature_dict.items():
            parent_masstrack_id = v['parent_masstrack_id']
            if parent_masstrack_id in mzstr_dict:
                mzstr_dict[parent_masstrack_id].append(n)
            else:
                mzstr_dict[parent_masstrack_id] = [n]

        return feature_dict, mzstr_dict

    def sort_nodes_by_mz(self):
        '''sort the nodes by increasing m/z
        '''
        return sorted( [(self.feature_dict[n]['mz'], n) for n in self.nodes_to_use ])
        
    def build_full_grid(self, isotopic_edges, adduct_edges):
        '''Build a khipu grid, after the input network is cleaned to isotopic_edges and adduct_edges.
        1) Get isotopic branches first, and treat them each branch as a node. This converts (U, V) to (U, B)
        2) Get minimum_spanning_tree, i.e. non-redundant adduct_edges from (U, B), 
            which should contain necessary adduct_edges.
        3) Order adducts as a vector, by topographic order. 
        4)  Generate grid starting from smallest m/z as (M0, A1).
        5) Match all nodes to the grid (snap_features_to_grid).

        Parameters
        ----------
        isotopic_edges, adduct_edges : separated from self.input_network

        Notes
        -----
        It's not a simple grid of combinations of isotope_search_patterns, adduct_search_patterns,
        because each pattern may occur more than once.
        Tree search methods don't always cover all edges, not a good choice.
        If minimum_spanning_tree is calculated from self.input_network, the isotopic branch may not be fully connect.

        Enforce unique node per grid, by best rtime with future option of a fitness function. 
        Extra nodes go to a new khipu.
        Designate feature of lowest m/z as root, which is often M+H+ or M-,  
        as we do not include lower mz in the initial search.
        A root is not necessarily M0, which may not be detected in perfect labeling experiments.
        '''
        abstracted_adduct_edges, root_branch, branch_dict = self.branch_abstraction(
            isotopic_edges, adduct_edges
        )
        indexed_adducts, adduct_index_labels, expected_grid_mz_values = self.build_grid_abstracted(
            abstracted_adduct_edges, root_branch
        )
        self.adduct_index = indexed_adducts
        self.adduct_index_labels = adduct_index_labels
        self.branch_dict = branch_dict
        self.khipu_grid = self.snap_features_to_grid(branch_dict, expected_grid_mz_values)

    def branch_abstraction(self, isotopic_edges, adduct_edges):
        '''Abstract a group of connecgted isotopic featrures into a branch.
        Reduce the input network to a set of abstracted_adduct_edges (hence new tree) of B-nodes.
        
        Returns
        -------
        abstracted_adduct_edges : list of nonredundant directed edges with data tag
        root_branch : ID of the branch where root feature (lowest m/z) resides
        branch_dict : dictionary, branch ID to member features/nodes.

        Notes
        -----
        Membership of B-nodes is returned as branch_dict, which is needed to realign to khipu grid.
        Without branch constraint, the grid realignment is error prone.
        '''
        G = nx.Graph(isotopic_edges)
        subnetworks = [G.subgraph(c).copy() for c in nx.connected_components(G)]
        _dict_branch, branch_dict = {}, {}
        _ii = 0
        for g in subnetworks:
            nodes = list(g.nodes())
            branch_id = "B" + str(_ii)
            branch_dict[branch_id] = nodes
            for n in nodes:
                _dict_branch[n] =  branch_id      # node membership is exclusive to subnetwork
            _ii += 1

        root_branch = _dict_branch.get(self.root, self.root)
        tmp, tracker, abstracted_adduct_edges = [], [], []
        for e in adduct_edges:
            if self.feature_dict[e[0]]['mz'] > self.feature_dict[e[1]]['mz']:
                e = (e[1], e[0], e[2])
            new_edge = (_dict_branch.get(e[0], e[0]), _dict_branch.get(e[1], e[1]), e[2])
            tmp.append(new_edge)

        # remove redundant edges
        for e in tmp:
            tt = make_edge_tag(e)
            if tt not in tracker:
                abstracted_adduct_edges.append(e)
                tracker.append(tt)

        return abstracted_adduct_edges, root_branch, branch_dict

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

    def build_branch_only_grid(self):
        '''When there's only a single adduct, this builds a branch for 'A1'.
        isotope_index is fixed based on inital isotope_search_patterns.
        isotopic_edges not needed as we re-align nodes.
        Updates self.khipu_grid
        '''
        _d = realign_isotopes(self.sorted_mz_peak_ids, self.isotope_search_patterns)
        self.khipu_grid = pd.DataFrame.from_dict(
            _d, orient="index", columns=['A1'], dtype=str,
        )
        for k,v in _d.items():
            self.annotation_dict[v] = (k, 'A1')
        
    def build_trunk_only_grid(self, adduct_edges):
        '''When no isotopes, only trunk is needed to describe adducts.
        Updates self.khipu_grid
        '''
        adduct_edges = nx.minimum_spanning_tree(nx.Graph(adduct_edges)).edges(data=True)
        dict_node_label = {}
        for e in adduct_edges:
            if self.feature_dict[e[0]]['mz'] > self.feature_dict[e[1]]['mz']:
                e = (e[1], e[0], e[2])
            dict_node_label[e[1]] = '(' + e[0] + '+' + e[2]['tag'] + ')'

        DG = nx.DiGraph(adduct_edges)
        indexed_adducts = list(nx.topological_sort(DG))     # ordered node IDs
        self.adduct_index = [dict_node_label.get(x, x) for x in indexed_adducts]
        self.khipu_grid = pd.DataFrame.from_dict(
            {"M": indexed_adducts}, 
            orient="index", columns=self.adduct_index, dtype=str,
        )
        for x in indexed_adducts:
            self.annotation_dict[x] = ('M', dict_node_label.get(x, x))

    def build_grid_abstracted(self, abstracted_adduct_edges, root_branch):
        '''Compute grid structure and expected_grid_mz_values.

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
        for e in nx.dfs_edges(DG, source=root_branch):
            indexed_adducts.append(e[1])
            target_mz = node_mz_dict[e[0]] + mz_modification_dict[make_edge_tag((e[0], e[1]))]
            # e[0] is guanranteed to be in node_mz_dict due to the tree walking order
            node_mz_dict[e[1]] = target_mz
            trunk_mzlist.append(target_mz)
            
        adduct_index_labels = [dict_node_label.get(x, x) for x in indexed_adducts]
        expected_grid_mz_values = []
        for A in trunk_mzlist:
            expected_grid_mz_values.append(
                [A] + [A+x[0] for x in self.isotope_search_patterns]
            ) 

        return indexed_adducts, adduct_index_labels, expected_grid_mz_values

    def extended_search(self, mztree, adduct_search_patterns_extended, mz_tolerance_ppm=5, rt_tolerance=2):
        '''Find additional adducts from unassigned_peaks using adduct_search_patterns_extended.
        mztree here are indexed unassigned_peaks.

        Updates
        -------
        self.annotation_dict
        self.khipu_grid
        self.adduct_index_labels

        Returns
        -------
        added_peaks

        Notes
        -----
        annotation_dict may have fewer nodes than nodes_to_use, but is preferred here for cleaner results.
        '''
        matched, _new_anno_dict, _grid_dict = [], {}, {}
        for n, v in self.annotation_dict.items():
            _iso, _ad = v
            P1 = self.feature_dict[n]
            for _pair in adduct_search_patterns_extended:
                (mass_difference, relation) = _pair[:2]
                _match = []
                tmp = find_all_matches_centurion_indexed_list(P1['mz'] + mass_difference, 
                                mztree, mz_tolerance_ppm)
                for P2 in tmp:
                    delta_rtime = abs(P1['rtime']-P2['rtime'])
                    if delta_rtime <= rt_tolerance:
                        _match.append((delta_rtime, P2))
                if _match:      # get best match by rtime only 
                    best_match_peak = sorted(_match)[0][1]
                    e1 = best_match_peak['id']
                    matched.append( (n, e1, relation) )
                    # P2 = P1 + adduct_relation. ',' avoids overwriting of existing names in khipu_grid
                    _new_anno_dict[e1] = (_iso, _ad + ',' + relation)
                    self.feature_dict[e1] = best_match_peak

        _new_adduct_index = list(set([x[1] for x in _new_anno_dict.values()]))
        _new_df = pd.DataFrame( np.empty((len(self.isotope_index), len(_new_adduct_index)), dtype=np.str),
                            index=self.isotope_index,
                            columns=_new_adduct_index,
                            dtype=str)
        for k,v in _new_anno_dict.items():
            _new_df.loc[v[0], v[1]] = k

        self.khipu_grid = pd.concat([self.khipu_grid, _new_df], axis=1)
        self.annotation_dict.update(_new_anno_dict)
        self.adduct_index_labels = list(self.khipu_grid.columns)

        return [x[1] for x in matched]

    def get_khipu_intensities(self):
        '''Return abundance_matrix as DataFrame in same layout as self.khipu_grid
        '''
        (_M, _N) = self.khipu_grid.shape
        abundance_matrix = pd.DataFrame(data=np.zeros(shape=(_M, _N)),
                            index=self.isotope_index,
                            columns=self.adduct_index_labels,
        )
        for ii in range(_M):
            for jj in range(_N):
                if self.khipu_grid.iloc[ii,jj]:
                    abundance_matrix.iloc[ii, jj] = self.feature_dict[self.khipu_grid.iloc[ii,jj]
                                                                ]['representative_intensity']
        return abundance_matrix


    #---------- export and visual --------------
    def print_khipu(self):
        '''Print khipu using adducts as trunk and isotopes as branches (preferred)
        '''
        print(self.khipu_grid)
        
    def print_khipu_rotated(self):
        '''Print khipu using isotopes as trunk and adducts as branches
        '''
        print(self.khipu_grid.T)
        
    def plot_khipu_diagram(self):
        '''Plot the khipu grid as diagram.
        Use MatPlotLib as default engine.
        '''
        df = self.get_khipu_intensities()
        _M, _N = df.shape
        zdata = []
        for ii in range(_M):
            for jj in range(_N):
                zdata.append((jj, ii, df.iloc[ii, jj]))
        
        X = [d[0] for d in zdata]
        Y = [d[1] for d in zdata]
        S = [(np.log10(d[2]+1))**2 for d in zdata]
        
        fig, ax = plt.subplots()
        
        for jj in range(_N):
            _t = str(int(df.iloc[0, jj])) + ' ~ ' + df.columns[jj]
            ax.text(jj, -1, _t, rotation=60)
            ax.plot([jj]*_M, range(_M), marker='o', linestyle='--', markersize=0.1)
        
        ax.plot([-1, _N+1], [0,0], linestyle='-', linewidth=2, color='k', alpha=0.3)
        ax.scatter(X, Y, c='red', s=S, alpha=0.8)
        
        for ii in range(_M):
            ax.text(_N+1.6, ii, df.index[ii])
        
        ax.margins(0.2)
        ax.set_axis_off()
        ax.invert_yaxis()
        
        #fig.tight_layout()
        plt.show()


    def plot_khipu_diagram_rotated(self):
        pass


    def export_json(self):
        '''Placeholder.
        '''
        return self.khipu_grid.to_json()


    def export_json_grid_data(self):
        return self.khipu_grid.to_json()


    def format_to_epds(self, id=''):
        '''Format khipu to empirical compound, with added ion notions.
        A small number of features do not map the khipu grid (e.g. their edges violate DAG rules). 
        They are kept in empCpd as "undetermined".
        '''
        if not id:
            id = 'root@' + self.root
        features = []
        for n in self.nodes_to_use:
            try:
                m, a = self.annotation_dict[n]
                ion = self.feature_dict[n]
                ion['isotope'] = m
                ion['modification'] = a
                ion['ion_relation'] = ','.join([m, a])
            except KeyError:
                # print("\nannotation_dict KeyError ", n, id)
                # print(self.clean_network.edges(data=True))
                # self.print()
                ion = self.feature_dict[n]
                ion['ion_relation'] = 'undetermined'
            features.append( ion )

        return {
            'interim_id': id, 
            'neutral_formula_mass': None,
            'neutral_formula': None,
            'Database_referred': [],
            'identity': [],
            'MS1_pseudo_Spectra': features,
            'MS2_Spectra': [],
            }


    #---------- shorthand names --------------
    print = print_khipu
    print2 = print_khipu_rotated
    plot = plot_khipu_diagram
    plot2 = plot_khipu_diagram_rotated

    #---------- placeholders --------------
    def select_by_fitness(self):
        '''After a khipu frame is built, among redundant features, 
        select the ones of best fitness score for the khipu.        
        '''
        pass

    def fitness(self):
        '''Fitness function of the khipu. 
        More a placeholder for now. Because unique assignment of a feature to a khipu can be based on
        a) closest retention time, and b) similar abundance patterns among adducts.
        The a) is depdendent on how well the samples were analyzed and data were preprocessed.
        The b) is not reliable as a pair of ions from another compound can still get good correlation 
        by disrupting the adduct patterns together.
        Default to a) is good enough for now.
        '''
        pass
