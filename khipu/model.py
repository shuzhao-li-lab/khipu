'''
Comprehensive construction of empCpds via generic tree structures.
Each annoTree = empCpd["MS1_pseudo_Spectra"].

        subnetwork : undirected graph. Example edges:
            [('F1606', 'F20', {'type': 'modification', 'tag': 'H'}), 
            ('F3533', 'F20', {'type': 'modification', 'tag': 'Na/H'}), 
            ('F195', 'F20', {'type': 'modification', 'tag': 'Acetonitrile'}), 
            ('F20', 'F807', {'type': 'modification', 'tag': 'Acetonitrile'}), 
            ('F20', 'F53', {'type': 'isotope', 'tag': '13C/12C'}), 
            ('F874', 'F808', {'type': 'isotope', 'tag': '13C/12C'})]

        edge_dict : edge_tag: {"source": e[0], "target": e[1], "type": "isotope", "tag": e[2]}
            edge_tag is str sorted, but the dict values preserve the direction, which is missed in nx.subnetwork.



    All isotopes of the same adduct form a khipu_string (branch). 
    All khipu_strings are attached to M0 khipu_rope (trunk).

        self.khipu_grid = {                 # JSON representation of khipu_grid
            'adduct_index': [],
            'isotope_index': [],
            'khipu_matrix': [],
        }


Trees use peak IDs, but full peak data are in json peak lists.


# Above only includes isotopic peaks; will add peaks without C13 counterparts into empCpd trees

# from itertools import combinations
# import networkx as nx
import treelib
'''

import json
import numpy as np
import pandas as pd

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
        '''
        self.input_network = subnetwork
        self.isotope_search_patterns = isotope_search_patterns
        self.adduct_search_patterns = adduct_search_patterns

        self.adduct_index = []
        self.isotope_index = ['M0'] + [x[1] for x in isotope_search_patterns]
        self.khipu_grid = []                # DataFrame of N x M, M number of isotopes, N number of adducts
                                            # M is explicitly specified during search
                                            # N is dynamically determined on each khipu

        self.pruned_network = None        
        self.is_13C_verified = False
        self.is_singleton = False

    def build_khipu(self, peak_dict, mz_tolerance_ppm=5, check_fitness=False):
        '''Convert a network of two types of edges, isotopic and adduct, to khipu instances of unique nodes.
        Use the grid to enforce unique ion in each position, as the initial network can contain 
        erroneously redundant leaves/nodes.

        check_fitness : if True, chech if replacing with any redundant nodes improves fitness. 
            This is a placeholder as a good fitness function is not easy to implement and it slows down computing.

        Updates
        -------
        self.khipu_grid : pd.DataFrame
        
        '''
        if self.input_network.size() == 2:
            self.build_simple_pair_grid(peak_dict)
        else:
            self.clean(peak_dict, mz_tolerance_ppm)
            isotopic_edges, adduct_edges = [], []
            for e in self.input_network.edges(data=True):
                if e[2]['type'] == 'isotope':
                    isotopic_edges.append(e)
                else:           # no other type of connections are allowed 
                    adduct_edges.append(e)
            
            if isotopic_edges and adduct_edges:
                self.build_full_grid(isotopic_edges, adduct_edges)
            elif isotopic_edges :       # branch only A1 adduct
                self.build_branch_only_grid(isotopic_edges)
            else:                       # trunk only, no isotopes
                self.build_trunk_only_grid(adduct_edges)
            
            if check_fitness:
                self.select_by_fitness()

            self.pruned_network = self.input_network[self.redundant_nodes]

    def build_simple_pair_grid(self, peak_dict):
        '''A khipu of only two nodes (one edge) does not need to go through full grid process,
        and is formatted here. The full grid process would work for these simple cases, but less efficient.
        '''
        e = self.input_network.edges()[0]
        if peak_dict[e[0]]['mz'] > peak_dict[e[1]]['mz']:
            e[0], e[1] = e[1], e[0]
        if e[2]['type'] == 'isotope':
            self.khipu_grid = pd.DataFrame( {'A1': [e[0], e[1]]},
                            index=['M0', e[2]['tag']], 
                            dtype=str)
        else:           # adduct
            self.khipu_grid = pd.DataFrame( {'A1': [e[0]], e[2]['tag']: [e[1]]},
                            index=['M0',], 
                            dtype=str)


    def clean(self, peak_dict, mz_tolerance_ppm):
        '''Clean up the input subnetwork, only using unique features to build a khipu frame.
        Redundant features are kept aside. 
        The leftover features are sent off to build new khipus.
        '''
        self.feature_dict, self.mzstr_dict = self.get_feature_dict(peak_dict, mz_tolerance_ppm)
        self.nodes_to_use = []
        self.redundant_nodes = []
        
        self.median_rtime = np.median([[(self.feature_dict[n]['rtime'], n) for n in self.input_network]])
        for k,v in self.mzstr_dict.items():
            if len(v) == 1:
                self.nodes_to_use.append(v[0]['id'])
            else:
                sorted_by_rtime_match = sorted(
                    [(abs(n['rtime'] - self.median_rtime), n['id']) for n in v]
                )
                self.nodes_to_use.append(sorted_by_rtime_match[0][1])
                self.redundant_nodes += [n[1] for n in sorted_by_rtime_match[1:]]

        self.sorted_mz_peak_ids = self.sort_nodes_by_mz()
        self.root = self.sorted_mz_peak_ids[0][1]


    def get_feature_dict(self, peak_dict, mz_tolerance_ppm):
        '''Index all input features; establish str identifier for features of same/close m/z values.
        Base on asari mass track IDs; keep unique m/z only.
        It's more efficient to use, since feature_dict is much smaller than peak_dict

        Parameters
        ----------
        peak_dict : dict of peaks/features indexed by IDs. Must have fields 
                    'id', 'mz', 'rtime', 'representative_intensity'
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
        return sorted( [(self.feature_dict[n]['mz'], n) for n in self.input_network.nodes() ])
        

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
        
        indexed_adducts, adduct_index_labels, expected_grid_mz_values = self.build_trunk_abstracted(
            self.branch_abstraction(isotopic_edges, adduct_edges)
        )
        self.adduct_index = adduct_index_labels
        self.khipu_grid = self.snap_features_to_grid(expected_grid_mz_values)


    def branch_abstraction(self, isotopic_edges, adduct_edges):
        '''Abstract a group of connecgted isotopic featrures into a branch.
        Reduce the input network to a set of abstracted_adduct_edges (hence new tree) of B-nodes.
        Membership of B-nodes is not returned, because all the feature nodes will be realigned to khipu grid.

        Returns
        -------
        abstracted_adduct_edges : list of nonredundant directed edges with data tag
        root_branch : ID of the branch where root feature (lowest m/z) resides
        '''
        G = nx.Graph(isotopic_edges)
        subnetworks = [G.subgraph(c).copy() for c in nx.connected_components(G)]
        _dict_branch = {}
        _ii = 0
        for g in subnetworks:
            for n in g.nodes():
                _dict_branch[n] = _ii       # node membership is exclusive to subnetwork
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

        return abstracted_adduct_edges, root_branch


    def snap_features_to_grid(self, expected_grid_mz_values):
        '''Create khipu_grid. To snap each feature to the expected_grid_mz_values.
        Parameters
        ----------
        expected_grid_mz_values : m/z values in a list of branches.

        Returns
        -------
        khipu_grid : pd.DataFrame, adducts as cols (trunk) and isotopes as rows (branches).

        Notes
        -----
        DataFrame is better for khipu_grid, because easier to cast data types and keep matrix format.
        Nested list or numpy array is harder to ensure correct matrix format.
        '''
        expected = np.array(expected_grid_mz_values)
        khipu_grid = pd.DataFrame( np.empty(expected.shape, dtype=np.str),
                            dtype=str)
        for x in self.sorted_mz_peak_ids:
            ii = np.argmin(abs(expected - x[0]))
            khipu_grid.iloc[np.unravel_index(ii, expected.shape)] = x[1]

        return khipu_grid


    def build_branch_only_grid(self, isotopic_edges):
        '''When there's only a single adduct, this builds a branch for 'A1'.
        isotope_index is fixed based on inital isotope_search_patterns.
        Updates
        -------
        khipu grid, including self.adduct_index, self.khipu_array
        '''
        nodes = set([])
        for e in isotopic_edges:
            nodes.add(e[0])
            nodes.add(e[1])
        sorted_mz_peak_ids = sorted([(self.feature_dict[n]['mz'], n) for n in nodes])
        _d = realign_isotopes(sorted_mz_peak_ids, self.isotope_search_patterns)
        self.khipu_grid = pd.DataFrame.from_dict(
            _d, orient="index", columns=['A1'], dtype=str,
        )
        

    def build_trunk_only_grid(self, adduct_edges):
        '''When no isotopes, only trunk is needed to describe adducts.
        '''
        adduct_edges = nx.minimum_spanning_tree(nx.Graph(adduct_edges)).edges(data=True)
        dict_node_label = {}
        for e in adduct_edges:
            if self.feature_dict[e[0]]['mz'] > self.feature_dict[e[1]]['mz']:
                e = (e[1], e[0], e[2])
            dict_node_label[e[1]] = e[1] + ' = ' + e[0] + ' + ' + e[2]['tag']

        DG = nx.DiGraph(adduct_edges)
        indexed_adducts = list(nx.topological_sort(DG))     # ordered node IDs
        self.adduct_index = [dict_node_label.get(x, x) for x in indexed_adducts]
        self.khipu_grid = pd.DataFrame.from_dict(
            {"M": indexed_adducts}, 
            orient="index", columns=self.adduct_index, dtype=str,
        )



    def build_grid_abstracted(self, abstracted_adduct_edges, root_branch):
        '''Compute grid structure and expected_grid_mz_values.

        Parameters
        ----------
        abstracted_adduct_edges : directed edges that connect B-nodes based on adduct m/z difference not isotopes.
            B-nodes are abstracted branches from isotopes during grid building.

        Returns
        -------
        indexed_adducts : a list containing ordered adducts
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

        DG = nx.DiGraph(abstracted_adduct_edges)
        # walk the graph through all nodes
        indexed_adducts, trunk_mzlist = [root_branch, ], [root_mz, ]
        for e in nx.dfs_edges(DG, source=root_branch):
            indexed_adducts.append(e[1])
            target_mz = node_mz_dict[e[0]] + mz_modification_dict[make_edge_tag((e[0], e[1]))]
            node_mz_dict[e[1]] = target_mz
            trunk_mzlist.append(target_mz)
            dict_node_label[e[1]] = e[1] + ' = ' + e[0] + ' + ' + e[2]['tag']

        adduct_index_labels = [dict_node_label.get(x, x) for x in indexed_adducts]
        expected_grid_mz_values = []
        for A in trunk_mzlist:
            expected_grid_mz_values.append(
                [A] + [A+x[0] for x in self.isotope_search_patterns]
            ) 

        return indexed_adducts, adduct_index_labels, expected_grid_mz_values


    def print_khipu(self):
        '''
        s = '\t' + '\t'.join(self.khipu_label_list) + '\n'
        for ii in range(len(self.isotope_index)):
            s += self.isotope_index + '\t' + '\t'.join(self.khipu_label_list) + '\n'
        '''
        print(self.khipu_grid)
        

    def plot_khipu_diagram(self):
        '''Plot the khipu grid as diagram.

        Use MatPlotLib as default engine?

        A khipu grid is a list of isotopic_branch dictionaries. 
        The list index has a companion list of Adduct labels.

        Use self.sorted_mz_peak_ids. Create node labels as Feature+Adduct.   ??

        '''
        


        self.abundance_matrix = np.zeros((len(self.isotope_index), len(self.adduct_index)))

        for f in self.feature_positions:
            self.abundance_matrix[f] = self.feature_dict[f]['representative intensity...            ']






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


    def export_json(self):
        return self.khipu_grid.to_json()


    def format_to_epds(self):
        '''
        Add to 'ion_relation' - 'isotope': [M0, C13*M3..], 'modification': ['-CO2', ...]


        '''
        pass


