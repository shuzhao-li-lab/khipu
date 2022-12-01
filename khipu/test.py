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


'''

import treelib
from mass2chem.io import read_features

from .model import khipu
from .utils import peaks_to_networks

adduct_search_patterns = [(1.0078, 'H'), (21.9820, 'Na/H'), (41.026549, 'Acetonitrile')]
isotope_search_patterns = [ (1.003355, '13C/12C', (0, 0.8)),
                                        (2.00671, '13C/12C*2', (0, 0.8)),
                                        (3.010065, '13C/12C*3', (0, 0.8)),
                                        (4.01342, '13C/12C*4', (0, 0.8)),
                                        (5.016775, '13C/12C*5', (0, 0.8)),
                                        (6.02013, '13C/12C*6', (0, 0.8)),]

def test_read(f='../testdata/full_Feature_table.tsv'):
    flist = read_features(f, id_col=0, mz_col=1, rtime_col=2, intensity_cols=(11, 17))
    subnetworks, peak_dict, edge_dict = peaks_to_networks(flist)
    return subnetworks, peak_dict, edge_dict



class khipu_interactive(khipu):

    def plot(self):
        pass


    def to_dataframe(self):
        pass


        
    def build_diagnostic_tree(self, peak_dict, depth_limit = 10):
        '''Build diagnostic tree before assigning nodes to grid and khipu diagram.

        Parameters
        ----------
        peak_dict : peaks (features) indexed by id.
        self.minimum_spanning_tree : minimum_spanning_tree from initial subnetwork

        Updates
        -------
        self.diagnostic_tree as an instance of treelib.Tree.

        Examples (tree here is diagnostic_tree)
        --------
        >>> from mass2chem.io import read_features
        >>> from khipu import model
        >>> f = 'testdata/full_Feature_table.tsv'
        >>> flist = read_features(f, id_col=0, mz_col=1, rtime_col=2, intensity_cols=(11, 17))
        table headers ordered:  mz rtime
        Read 4016 feature lines
        >>> subnetworks, peak_dict, edge_dict = model.peaks_to_networks(flist,)
        >>> big = [g for g in subnetworks if g.size()>5]
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
        tree = treelib.Tree()
        root_node = self.root


        T = self.minimum_spanning_tree


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

        self.diagnostic_tree = tree
        if remaining:
            print(remaining)


