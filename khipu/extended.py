import treelib
from .model import *

class khipu_diagnosis(khipu):
    '''Added diagnostic functions to khipu class.
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
        root_node = int(peak_dict[self.root]['mz']) - 1
        tree.create_node(str(root_node), root_node)
        # tree depth = 2, defined by khipu
        for n in self.adduct_index:
            pdata = peak_dict.get(n, {})
            plabel = n
            if pdata:
                plabel = make_peak_tag(pdata)
            tree.create_node(plabel, n, parent=root_node, data=pdata)
            if n in self.branch_dict:
                for b in self.branch_dict[n]:
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

    def plot(self):
        pass

    def to_dataframe(self):
        pass
