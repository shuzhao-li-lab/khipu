'''
>>> from mass2chem.io import read_features
>>> 
>>> f = 'testdata/full_Feature_table.tsv'
>>> 
>>> flist = read_features(f, id_col=0, mz_col=1, rtime_col=2, intensity_cols=(11, 17))
table headers ordered:  mz rtime
Read 4016 feature lines
>>> 
>>> print(len(flist), flist[5])
4016 {'id': 'F6', 'mz': 117.0659, 'rtime': 148.54, 'intensities': [5630973.0, 321067.0, 237998.0, 97759.0], 'representative_intensity': 1571949.25}
>>> 

>>> from khipu.test import *

>>> subnetworks, peak_dict, edge_dict = test_read('testdata/full_Feature_table.tsv')
table headers ordered:  mz rtime
Read 4016 feature lines
>>> 

>>> subnetworks, peak_dict, edge_dict = test_read_url()
Retrieving test data from GitHub.
table headers ordered:  mz rtime
Read 4016 feature lines

>>> big = [g for g in subnetworks if g.size()>8]
>>> 
>>> KP = khipu(big[1], isotope_search_patterns, adduct_search_patterns)
>>> 
>>> KP.build_khipu(peak_dict)
>>> 
>>> KP.print_khipu()
              B0 B1 = F1606 + H
M0          F195            F20
13C/12C                     F53
13C/12C*2                      
13C/12C*3                      
13C/12C*4                  F254
13C/12C*5                 F3271
13C/12C*6  F3398           F874


from pathlib import Path
import sys
path_root = Path(__file__).parents[2]
sys.path.append(str(path_root))
print(sys.path)

➜  khipu git:(main) ✗ python3 -m khipu.test


'''
import urllib.request
import treelib
import numpy as np

from .model import khipu
from .utils import *

adduct_search_patterns = [(1.0078, 'H'), (21.9820, 'Na/H'), (41.026549, 'Acetonitrile')]
isotope_search_patterns = [ (1.003355, '13C/12C', (0, 0.8)),
                                        (2.00671, '13C/12C*2', (0, 0.8)),
                                        (3.010065, '13C/12C*3', (0, 0.8)),
                                        (4.01342, '13C/12C*4', (0, 0.8)),
                                        (5.016775, '13C/12C*5', (0, 0.8)),
                                        (6.02013, '13C/12C*6', (0, 0.8)),]


def test_read_file(f='../testdata/full_Feature_table.tsv'):
    '''Use from mass2chem.io import read_features
    '''
    flist = read_features(f, id_col=0, mz_col=1, rtime_col=2, intensity_cols=(11, 17))
    subnetworks, peak_dict, edge_dict = peaks_to_networks(flist)
    return subnetworks, peak_dict, edge_dict

def test_read_url(url='https://github.com/shuzhao-li/khipu/raw/main/testdata/full_Feature_table.tsv'):
    print("Retrieving test data from GitHub.")
    data = urllib.request.urlopen(url).read().decode('utf-8')
    flist = read_features_from_text(data, id_col=0, mz_col=1, rtime_col=2, intensity_cols=(11, 17))
    subnetworks, peak_dict, edge_dict = peaks_to_networks(flist)
    return subnetworks, peak_dict, edge_dict

def read_features_from_text(text_table, 
                        id_col=0, mz_col=1, rtime_col=2, 
                        intensity_cols=(3,4), delimiter="\t"):
    '''
    Read a text feature table into a list of features.
    Input
    -----
    text_table: Tab delimited feature table read as text. First line as header.
                    Recommended col 0 for ID, col 1 for m/z, col 2 for rtime.
    id_col: column for id. If feature ID is not given, row_number is used as ID.
    mz_col: column for m/z.
    rtime_col: column for retention time.
    intensity_cols: range of columns for intensity values. E.g. (3,5) includes only col 3 and 4.
    Return
    ------
    List of features: [{'id': '', 'mz': 0, 'rtime': 0, 
                        intensities: [], 'representative_intensity': 0, ...}, 
                        ...], 
                        where representative_intensity is mean value.
    '''
    # featureLines = open(feature_table).read().splitlines()
    featureLines = text_table.splitlines()
    header = featureLines[0].split(delimiter)
    num_features = len(featureLines)-1
    # sanity check
    print("table headers ordered: ", header[mz_col], header[rtime_col])
    print("Read %d feature lines" %num_features)
    L = []
    for ii in range(1, num_features+1):
        if featureLines[ii].strip():
            a = featureLines[ii].split(delimiter)
            if isinstance(id_col, int):         # feature id specified
                iid = a[id_col]
            else:
                iid = 'row'+str(ii)
            xstart, xend = intensity_cols
            intensities = [float(x) for x in a[xstart: xend]]
            L.append({
                'id': iid, 'mz': float(a[mz_col]), 'rtime': float(a[rtime_col]),
                'intensities': intensities,
                'representative_intensity': np.mean(intensities),
            })
    return L


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


#
# -----------------------------------------------------------------------------
#
if __name__ == '__main__':
    subnetworks, peak_dict, edge_dict = test_read_url()
    big = [g for g in subnetworks if g.size()>8]

    print("\n\n")
    print("Example khipu of one empirical compound from demo data.")
    KP = khipu_diagnosis(big[ -4 ], isotope_search_patterns, adduct_search_patterns)
    print(KP.input_network.edges(data=True))
    KP.build_khipu(peak_dict)
    print(KP.sorted_mz_peak_ids, "\n")
    KP.show_trimming()
    print("\n\n")
    KP.build_diagnostic_tree_full()

    print("Looking at khipu: ")
    KP.print()
    print("\n\n")
    print("Rotated way to look at khipu: ")
    KP.print2()
    print("\n\n")
    KP.build_diagnostic_tree_clean()

    KP.plot_khipu_diagram()

    print("\n\n")
    print("Multiple example khipus: ")
    print("======================== \n")
    for g in big[ -10: ]:
        KP = khipu(g, isotope_search_patterns, adduct_search_patterns)
        KP.build_khipu(peak_dict)
        print(KP.sorted_mz_peak_ids, "\n")
        KP.print_khipu()
        print('\n\n')

