'''
To run test from top directory:
➜ python3 -m khipu.test
'''
import urllib.request
import numpy as np

from .model import khipu
from .extended import *
from .utils import *


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
