'''
Comprehensive construction of empCpds via generic tree structures (khipu).
Each khipu = empCpd["MS1_pseudo_Spectra"].
This module contains two classes:
Weavor is constructed once and provides the main algorithms and functions;
Khipu is mainlty the data structures and built as one instance per emprical compound.
So from one experiment/feature table, one may get 10**3 khipus.
'''

import numpy as np
import pandas as pd
import treelib
from scipy.optimize import curve_fit
from .plot import plot_khipugram
from .utils import *


class Weavor:
    '''For each experiment, this class sets up a grid of isotopes and adducts,
    and provide function solve to produce matched features and inferred neutral mass.
    '''
    def __init__(self, peak_dict, 
                isotope_search_patterns, adduct_search_patterns, 
                mz_tolerance_ppm=5, mode='pos'):
        '''Initiate mzgrid and indices.        
        '''
        self.mode = mode
        self.mz_tolerance_ppm = mz_tolerance_ppm
        self.peak_dict = peak_dict
        self.isotope_search_patterns  = sorted(isotope_search_patterns)   # m/z sort search patterns
        self.adduct_search_patterns = sorted(adduct_search_patterns)
        self.make_grid()

    def make_grid(self):
        '''Create a grid of m/z values as DataFrame.
        adduct_pattern is computed using neutral mass as 0 offset.
        The orders of isotope_search_patterns and adduct_search_patterns are kept in grid,
        to enable easy mapping of edges to the grid.

        Updates
        -------
        self.isotope_index
        self.isotope_dict
        self.adduct_pattern
        self.adduct_dict
        self._size, _M, _N
        self.mzgrid : pd.DataFrame with isotopes as rows and adducts as cols
        '''
        _d, self.adduct_dict, self.isotope_dict = {}, {}, {}
        self.isotope_index = ['M0'] + [x[1] for x in self.isotope_search_patterns]
        for x in self.isotope_search_patterns:
            self.isotope_dict[x[1]] = x[0]

        self.adduct_pattern = make_expected_adduct_index(self.mode, self.adduct_search_patterns)
        # [(PRONTON, M+H+), ..., (42.0338, ACN), ...]
        self.adduct_index = [A[1] for A in self.adduct_pattern]
        for A in self.adduct_pattern:
            _d[A[1]] = [A[0]] + [A[0]+x[0] for x in self.isotope_search_patterns]
            self.adduct_dict[A[1]] = A[0]

        self._M, self._N = len(self.isotope_index), len(self.adduct_index)
        self._size = self._M * self._N 
        self.mzgrid = pd.DataFrame(_d, index=self.isotope_index)       

    def build_simple_pair_grid(self, e):
        '''A khipu of only two nodes (one edge) does not need to go through full grid process,
        and is formatted here. The full grid process would work for these simple cases, but less efficient.
        neutral_mass is taken as average of two fitted values.
        '''
        mz1, mz2 = self.peak_dict[e[0]]['mz'], self.peak_dict[e[1]]['mz']
        if mz1 > mz2:
            e = (e[1], e[0], e[2])
            mz1, mz2 = mz2, mz1

        if e[2]['type'] == 'isotope':
            grid = pd.DataFrame( {self.adduct_index[0]: [e[0], e[1]]},
                            index=[self.isotope_index[0], e[2]['tag']], 
                            dtype=str)
            feature_map = {
                e[0]: (self.isotope_index[0], self.adduct_index[0]), 
                e[1]: (e[2]['tag'], self.adduct_index[0]),
            }
            
        else:           # adduct
            grid = pd.DataFrame( {self.adduct_index[0]: [e[0]], e[2]['tag']: [e[1]]},
                            index=[self.isotope_index[0],], 
                            dtype=str)
            feature_map = {
                e[0]: (self.isotope_index[0], self.adduct_index[0]), 
                e[1]: (self.isotope_index[0], e[2]['tag']),
            }

        neutral_mass = self.regress_neutral_mass(feature_map)
        return neutral_mass, grid, feature_map

    def build_branch_only_grid(self, sorted_mz_peak_ids):
        '''When there's only a single adduct, this builds a branch for adduct_index[0].
        isotope_index is fixed based on inital isotope_search_patterns.
        returns neutral_mass, grid
        '''
        _d = realign_isotopes(sorted_mz_peak_ids, self.isotope_search_patterns)
        # {'M0': F1, '13C/12C*2': F11, ...}
        grid = pd.DataFrame.from_dict(
            _d, orient="index", columns=[self.adduct_index[0]], dtype=str,
        )
        feature_map = {}
        for k,v in _d.items():
            feature_map[v] = (k, self.adduct_index[0])
        neutral_mass = self.regress_neutral_mass(feature_map)
        #neutral_mass = sorted_mz_peak_ids[0][0] - self.adduct_pattern[0][0] # relative to M+H+ or M-H-

        return neutral_mass, grid, feature_map
        
    def build_trunk_only_grid(self, adduct_edges):
        '''When no isotopes, only trunk is needed to describe adducts.
        returns neutral_mass, grid
        '''
        root, edges = self.trunk_solver(adduct_edges)
        _d = {self.adduct_index[0]: root}
        for e in edges:
            _d[e[2]['tag']] = e[1]
        grid = pd.DataFrame(_d, index=[self.isotope_index[0],])
        feature_map = {}
        for k,v in _d.items():
            feature_map[v] = (self.isotope_index[0], k)
        neutral_mass = self.regress_neutral_mass(feature_map)

        return neutral_mass, grid, feature_map

    def regress_neutral_mass(self, feature_map):
        '''Get neutral mass by regression on mapped features.
        feature_map : {feature_id: (isotope_index, adduct_index)}
        return 
        '''
        def _func(x, neu):
            return x + neu
        fixed_feature_list = list(feature_map.keys())
        Y = [self.peak_dict[f]['mz'] for f in fixed_feature_list]
        X = [self.mzgrid.at[feature_map[f][0], feature_map[f][1]] for f in fixed_feature_list]
        popt, _ = curve_fit(_func, X, Y)
        neutral_mass = popt[0]
        return neutral_mass

    def trunk_solver(self, adduct_edges, branch_dict={}):
        '''Find best solution to fit a set of edges on the trunk. 
        This uses score_graph_on_trunk to score matched graph, which can be 
        abstracted_adduct_edges or regular adduct edges.

        adduct_edges : m/z ordered edge with tag on edge ion relationship, [(n1, n2, relation), ...].
        branch_dict : when abstract branch is used, this is the dict for memberships.

        returns selected best root and subset of edges.
        '''
        nodes = set([])
        for e in adduct_edges:
            nodes.add(e[0])
            nodes.add(e[1])
        weights = {}
        for x in nodes:
            if x in branch_dict:
                weights[x] = len(branch_dict[x])
            else:
                weights[x] = 1
        
        scored = sorted(
            [self.score_graph_on_trunk(root, adduct_edges, weights) for root in nodes]
        )
        _, root, edges = scored[-1]

        return root, edges

    def score_graph_on_trunk(self, root, adduct_edges, weights):
        '''Use weighted algorithm on reference adduct tree to calcualte match score.
        returns score, root, selected_edges
        '''
        selected_edges = []
        score = weights[root]
        for e in adduct_edges:
            if e[0] == root:
                selected_edges.append(e)
                score += weights[e[1]]

        return score, root, selected_edges

    def score_graph_on_grid(self, root_corrected_mz_features, mz_error=0.01):
        '''Check how many values in root_corrected_mzs match to self.mzgrid; count = score.
        returns score, feature_map
        '''
        feature_map = {}
        for mz, f in root_corrected_mz_features:
            delta = abs(self.mzgrid - mz)
            if delta.values.min() < mz_error:
                ii, jj = np.unravel_index( np.argmin(delta.values), self.mzgrid.shape )
                feature_map[f] = (self.mzgrid.index[ii], self.mzgrid.columns[jj])

        score = len(feature_map)

        return score, feature_map

    def build_full_grid(self, abstracted_adduct_edges, branch_dict, nodes_to_use):
        '''Build a khipu grid, after the input network is cleaned to isotopic_edges and adduct_edges.
        1) Get isotopic branches first, and treat them each branch as a node. This converts (U, V) to (U, B).
            Done in Khipu.branch_abstraction().
        2) Get optimal order of abstracted adduct_edges. Done in Weavor.trunk_solver(), by topology not involving m/z. 
        3) Generate grids starting from the smallest m/z in each branch. Get the grid of best feature map.
        4) Use the optimal feature map to set grids.

        Parameters
        ----------
        abstracted_adduct_edges : list of nonredundant directed edges with data tag.
            A node here can be a feature or a branch, which is a list of isotopes.
        branch_dict : dictionary, branch ID to member features/nodes.

        Notes
        -----
        We don't know the real root to start with, and can't assume the lowest m/z is M+H+ or M-H-.
        The root should be inferred from best overall pattern match.
        A pseudo-root is not necessarily M0, which may not be detected in perfect labeling experiments.

        Enforce unique node per grid, by best rtime with future option of a fitness function (Khipu.clean()). 
        Extra nodes go to a new khipu.
        '''
        mz_features = [(self.peak_dict[f]['mz'], f) for f in nodes_to_use]
        root, edges = self.trunk_solver(abstracted_adduct_edges, branch_dict)
        abstracted_adducts = [root, ]
        abstracted_adduct_tags = [self.adduct_index[0], ]
        for e in edges:
            abstracted_adducts.append(e[1])
            abstracted_adduct_tags.append(e[2]['tag'])

        list_bagged_features = []
        for x in abstracted_adducts:
            list_bagged_features.append(branch_dict.get(x, [x]))

        # this generates a list of possible neutral mass values
        peudo_roots = [0, ]         # 0 will not be used, just spacer
        for ii in range(len(abstracted_adducts)):
            L = list_bagged_features[ii]
            _offset = self.adduct_dict[abstracted_adduct_tags[ii]]
            peudo_roots.append(
                min([self.peak_dict[f]['mz'] for f in L]) - _offset
            )

        # will check which neutral mass is best fit
        list_grid_fits = []
        for ii in range(1, len(peudo_roots)):
            if abs(peudo_roots[ii] - peudo_roots[ii-1]) > 0.000001 * self.mz_tolerance_ppm * peudo_roots[ii]:
                score, feature_map = self.score_graph_on_grid(
                    [(x[0] - peudo_roots[ii], x[1]) for x in mz_features]
                )
                list_grid_fits.append((score, ii, feature_map))  # ii used as tie breaker in sorting

        best_feature_map = sorted(list_grid_fits, reverse=True)[0][2]
        neutral_mass = self.regress_neutral_mass(best_feature_map)

        grid = pd.DataFrame( np.empty((self._M, self._N), dtype=str),
                            index=self.mzgrid.index,
                            columns=self.mzgrid.columns,
                            dtype=str)
        for f,v in best_feature_map.items():
            grid.loc[v] = f

        return neutral_mass, grid, best_feature_map

    def make_tree(self):
        tree = treelib.Tree()
        tree.create_node('M', 'M')
        for node in [x[1] for x in self.adduct_pattern]:
            tree.create_node(node, node, parent='M')
            for m in self.isotope_index:
                tree.create_node(m, m+node, parent=node)

        return tree 

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


class Khipu:
    '''2-tier tree representation of an empirical compound.
    An empTree follows the matrix notion of a family of peaks: a tree of 2 levels - isotopes and adducts.
    A true root of khipu is the compound in neutral formula.
    A peudo root of khipu is the ion of lowest m/z.
    The representative adduct should have the most abundance.
    The matrix notion of a khipu is a DataFrame: columns are Adduct labels; rows isotopes.
    trunk : default for adducts
    branch : default for isotopes
    '''
    def __init__(self, subnetwork):
        '''Do not include negative m/z difference in this class, as the node of lowest m/z is used as root.
        Neutral loss and fragments can be searched after khipu is established.
        Input subnetwork is undirected graph, but khipu will require directed graph/tree, by increasing m/z.
        subnetwork : undirected graph. Example edges:
            ('F3533', 'F20', {'type': 'modification', 'tag': 'Na/H'}), 
            ('F195', 'F20', {'type': 'modification', 'tag': 'Acetonitrile'}), 
            ('F20', 'F53', {'type': 'isotope', 'tag': '13C/12C'}), 
        '''
        self.id = ''
        self.input_network = subnetwork
        self.nodes_to_use = []
        self.redundant_nodes = []          # nodes in input_network but not in final khipu
        self.pruned_network = None
        self.exceptions = []
        self.khipu_grid = {}                # DataFrame of N x M, M number of isotopes, N number of adducts
                                            # M is explicitly specified during search
                                            # N is dynamically determined on each khipu

    def build_khipu(self, WeavorInstance, mz_tolerance_ppm=5, check_fitness=False):
        '''Convert a network of two types of edges, isotopic and adduct, to khipu instances of unique nodes.
        Use the grid to enforce unique ion in each position, as the initial network can contain 
        erroneously redundant leaves/nodes.

        WeavorInstance : Weavor instance to solve the grid.
        check_fitness : if True, chech if replacing with any redundant nodes improves fitness. 
            This is a placeholder as a good fitness function is not easy to implement and it slows down computing.
            Future use can be if check_fitness: self.select_by_fitness()

        Updates
        -------
        self.khipu_grid : pd.DataFrame
        self.neutral_mass : inferred neutral mass for the khipu compound
        '''
        self._size_limit_ = WeavorInstance._size    # Set max limit of feature number based on grid size
        self.feature_dict, self.mzstr_dict = self.get_feature_dict(
                                        WeavorInstance.peak_dict, mz_tolerance_ppm)
        if self.input_network.number_of_edges() == 1:
            edge = list(self.input_network.edges(data=True))[0]
            self.neutral_mass, self.khipu_grid, self.feature_map =\
                        WeavorInstance.build_simple_pair_grid(edge)
            self.clean_network = self.input_network
            self.nodes_to_use = list(self.clean_network.nodes())

        else:
            self.clean(WeavorInstance, mz_tolerance_ppm)
            self.clean_network = self.input_network.subgraph(self.nodes_to_use)
            
            isotopic_edges, adduct_edges = [], []
            for e in self.clean_network.edges(data=True):
                if e[2]['type'] == 'isotope':
                    isotopic_edges.append(e)
                else:                   # adducts, no other type of edges are allowed 
                    adduct_edges.append(e)

            if isotopic_edges and adduct_edges:
                abstracted_adduct_edges, branch_dict = self.branch_abstraction(
                        isotopic_edges, adduct_edges
                        )  
                self.neutral_mass, self.khipu_grid, self.feature_map = WeavorInstance.build_full_grid(
                        abstracted_adduct_edges, branch_dict, self.nodes_to_use
                        )
            elif isotopic_edges :       # branch only single adduct
                self.neutral_mass, self.khipu_grid, self.feature_map = \
                    WeavorInstance.build_branch_only_grid(self.sorted_mz_peak_ids)
            elif adduct_edges:          # trunk only, no isotopes
                self.neutral_mass, self.khipu_grid, self.feature_map = \
                    WeavorInstance.build_trunk_only_grid(adduct_edges)
            else:
                print("Empty network - ", self.nodes_to_use, isotopic_edges, adduct_edges, 
                                        self.input_network.edges(data=True))
            
    def get_pruned_network(self):
        '''Get extra features and edges that are not fit in this khipu.
        Updates:
        self.redundant_nodes
        self.pruned_network
        '''
        self.redundant_nodes = [n for n in self.nodes_to_use if n not in self.feature_map]
        if self.redundant_nodes:
            self.pruned_network = self.input_network.subgraph(self.redundant_nodes)

    def clean(self, WeavorInstance, mz_tolerance_ppm):
        '''Clean up the input subnetwork, only using unique features to build a khipu frame.
        Redundant features are kept aside. The leftover features will be sent off to build new khipus.

        Updates
        -------
        self.feature_dict
        self.mzstr_dict
        self.median_rtime
        self.nodes_to_use
        self.redundant_nodes 
        self.sorted_mz_peak_ids
        self.root       # temporary

        Notes
        -----
        Note: mzstr_dict can be problematic in data that were not processed well, 
        because minor potential shift can confuse ion relations.
        This clean() may cut some initial edges.
        '''
        if self.input_network.number_of_nodes() > self._size_limit_:
            self.input_network = self.down_size()
        self.feature_dict, self.mzstr_dict = self.get_feature_dict(WeavorInstance.peak_dict, mz_tolerance_ppm)
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

    def down_size(self):
        '''If input_network is too big, down size by selected the features of highest abundance.
        The extra nodes are not kept, because they can be picked up by extended search later if they fit this Khipu.
        '''
        use_nodes = [
            (self.feature_dict[n]['representative_intensity'], n) for n in self.input_network.nodes()
        ]
        _N = len(use_nodes)
        use_nodes = sorted(use_nodes, reverse=True)[:self._size_limit_]
        print("Downsized input network with %d features, highest peak at %s " %(_N, use_nodes[0][1]))
        new_network = self.input_network.subgraph([x[1] for x in use_nodes])
        return new_network

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
        
    def branch_abstraction(self, isotopic_edges, adduct_edges):
        '''Abstract a group of connecgted isotopic featrures into a branch.
        Reduce the input network to a set of abstracted_adduct_edges (hence new tree) of B-nodes.
        
        Returns
        -------
        abstracted_adduct_edges : list of nonredundant directed edges with data tag
        branch_dict : dictionary, branch ID to member features/nodes.

        Notes
        -----
        Membership of B-nodes is returned as branch_dict, which is needed to realign to khipu grid.
        Without branch constraint, the grid realignment is error prone.
        Not checking if abstracted_adduct_edges are fully connected.
        '''
        G = nx.Graph(isotopic_edges)
        subnetworks = [G.subgraph(c).copy() for c in nx.connected_components(G)]
        _dict_branch, branch_dict = {}, {}
        _ii = 0
        for g in subnetworks:
            nodes = list(g.nodes())
            branch_id = "B" + str(_ii)            # id for internal use
            branch_dict[branch_id] = nodes
            for n in nodes:
                _dict_branch[n] =  branch_id      # node membership is exclusive to subnetwork
            _ii += 1

        # root_branch = _dict_branch.get(self.root, self.root)
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

        return abstracted_adduct_edges, branch_dict

    def extended_search(self, 
                        mztree, adduct_search_patterns_extended, 
                        mz_tolerance_ppm=5, rt_tolerance=2):
        '''Find additional adducts from unassigned_peaks using adduct_search_patterns_extended.
        mztree here are indexed unassigned_peaks.

        Updates
        -------
        self.feature_map
        self.khipu_grid

        Returns
        -------
        added_peaks

        Notes
        -----
        annotation_dict may have fewer nodes than nodes_to_use, but is preferred here for cleaner results.
        '''
        matched, _new_anno_dict = [], {}
        for n, v in self.feature_map.items():
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
                        _match.append((delta_rtime, P2['id'], P2))    # need id to break tie in sorting
                if _match:      # get best match by rtime only 
                    best_match_peak = sorted(_match)[0][-1]
                    e1 = best_match_peak['id']
                    matched.append( (n, e1, relation) )
                    # P2 = P1 + adduct_relation. ',' avoids overwriting of existing names in khipu_grid
                    _new_anno_dict[e1] = (_iso, _ad + ',' + relation)
                    self.feature_dict[e1] = best_match_peak

        _new_adduct_index = list(set([x[1] for x in _new_anno_dict.values()]))
        _new_df = pd.DataFrame( np.empty((self.khipu_grid.shape[0], len(_new_adduct_index)), 
                            dtype=str),
                            index=self.khipu_grid.index,
                            columns=_new_adduct_index,
                            dtype=str)
        for k,v in _new_anno_dict.items():
            _new_df.loc[v[0], v[1]] = k

        self.khipu_grid = pd.concat([self.khipu_grid, _new_df], axis=1)
        self.feature_map.update(_new_anno_dict)

        return [x[1] for x in matched]

    def get_khipu_intensities(self):
        '''Return abundance_matrix as DataFrame in same layout as self.khipu_grid
        '''
        (_M, _N) = self.khipu_grid.shape
        abundance_matrix = pd.DataFrame(data=np.zeros(shape=(_M, _N)),
                            index=self.khipu_grid.index,
                            columns=self.khipu_grid.columns,
        )
        for ii in range(_M):
            for jj in range(_N):
                if self.khipu_grid.iloc[ii,jj]:
                    abundance_matrix.iloc[ii, jj] = self.feature_dict[self.khipu_grid.iloc[ii,jj]
                                                                ]['representative_intensity']
        return abundance_matrix

    def get_khipu_mzgrid_print(self):
        '''Return str m/z matrix as DataFrame in same layout as self.khipu_grid, for visual purpose.
        '''
        _new_df = pd.DataFrame( np.empty(self.khipu_grid.shape, 
                            dtype=str),
                            index=self.khipu_grid.index,
                            columns=self.khipu_grid.columns,
                            dtype=str)
        for k,v in self.feature_map.items():
            _new_df.loc[v[0], v[1]] = self.feature_dict[k]['mz']

        return _new_df

    #---------- export and visual --------------
    def print_khipu(self):
        '''Print khipu using adducts as trunk and isotopes as branches (preferred)
        '''
        print(self.khipu_grid)
        
    def print_khipu_rotated(self):
        '''Print khipu using isotopes as trunk and adducts as branches
        '''
        print(self.khipu_grid.T)
        
    def plot_khipu_diagram(self, savepdf=''):
        '''Plot the khipu grid as diagram.
        Use MatPlotLib as default engine.
        '''
        plot_khipugram(self.get_khipu_intensities(), savepdf)

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
            id = 'root@' + str(round(self.neutral_mass, 4))
        features = []
        for n in self.nodes_to_use:
            try:
                m, a = self.feature_map[n]
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
            'neutral_formula_mass': self.neutral_mass,
            'neutral_formula': None,
            'Database_referred': [],
            'identity': [],
            'MS1_pseudo_Spectra': features,
            'MS2_Spectra': [],
            }


    #---------- shorthand names --------------
    print = print_khipu
    print2 = print_khipu_rotated
    print3 = get_khipu_mzgrid_print         # dataframe
    plot = plot_khipu_diagram
    plot2 = plot_khipu_diagram_rotated
    
