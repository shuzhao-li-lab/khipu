'''
To run test from top directory:
➜ python3 -m khipu.test
'''
import urllib.request
import numpy as np

from mass2chem.io import read_features

# from .json_encoder import NpEncoder
from .model import khipu
from .extended import *
from .utils import *


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
    subnetworks, peak_dict, edge_dict = peaks_to_networks(flist,
                        adduct_search_patterns = adduct_search_patterns,
                        isotope_search_patterns = isotope_search_patterns,
                        rt_tolerance = 1,       # 1 sec
    )
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


#
# -----------------------------------------------------------------------------
#
if __name__ == '__main__':
    subnetworks, peak_dict, edge_dict = test_read_url()
    big = [g for g in subnetworks if g.size()>15]

    print("\n\n")
    print("Example khipu of one empirical compound from demo data.")
    KP = khipu_diagnosis(big[ 1 ], isotope_search_patterns, adduct_search_patterns)
    print("Input network edges: ")
    print(KP.input_network.edges(data=True))
    KP.build_khipu(peak_dict)
    print(KP.sorted_mz_peak_ids, "\n")
    KP.show_trimming()
    print("\n\n")
    KP.build_diagnostic_tree_full()
    print("\n\n")
    KP.build_diagnostic_tree_clean()

    print("Looking at khipu: ")
    KP.print()
    print("\n\n")
    print("Rotated way to look at khipu: ")
    KP.print2()
    print("\n\n")
    KP.build_khipu_tree()

    unassigned_peaks = [v for x,v in peak_dict.items() if x not in KP.nodes_to_use]
    mztree =  build_centurion_tree(unassigned_peaks)
    KP.extended_search(mztree, extended_adducts)
    print("extended khipu with additional adducts: ")
    KP.print()

    # KP.plot_khipu_diagram()

    print("\n\n")
    print("Multiple example khipus: ")
    print("======================== \n")
    for g in big[ -6: ]:
        KP = khipu(g, isotope_search_patterns, adduct_search_patterns)
        KP.build_khipu(peak_dict)
        print(KP.sorted_mz_peak_ids, "\n")
        KP.print_khipu()
        print('\n\n')

    print("\n\n")
    print("Build and export JSON file of full list of khipus: ")
    print("================================================== \n")

    khipu_list = peak_dict_to_khipu_list(
        subnetworks, peak_dict, isotope_search_patterns, adduct_search_patterns)
    
    # khipu_list = extend_khipu_list(khipu_list, peak_dict, adduct_search_patterns_extended)
    print("\n\n ~~~~~~ Got %d khipus ~~~~~~~ \n\n" %len(khipu_list))
    empCpds = export_empCpd_khipu_list(khipu_list)

    outfile = 'khipu_test_empricalCompounds.json'
    with open(outfile, 'w', encoding='utf-8') as f:
        json.dump(empCpds, f, ensure_ascii=False, indent=2)
