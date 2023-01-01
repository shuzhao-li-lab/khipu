'''
To run test from top directory:
➜ python3 -m khipu.test
'''
import urllib.request
# import numpy as np

from khipu import __version__

from .model import Khipu
from .extended import *
from .utils import *

from .epdsConstructor import epdsConstructor

url= 'https://github.com/shuzhao-li/khipu/raw/main/testdata/ecoli_pos.tsv' 

def test_read_url(url=url):
    flist = get_peaklist_from_url(url)
    subnetworks, peak_dict, edge_dict = peaks_to_networks(flist,
                        adduct_search_patterns = adduct_search_patterns,
                        isotope_search_patterns = isotope_search_patterns,
                        rt_tolerance = 1,       # 1 sec
    )
    return subnetworks, peak_dict, edge_dict

def get_peaklist_from_url(url):
    print("Retrieving test data from GitHub.")
    data = urllib.request.urlopen(url).read().decode('utf-8')
    flist = read_features_from_text(data, id_col=0, mz_col=1, rtime_col=2, intensity_cols=(3, 9))
    return flist

def test1():
    '''Test khipu construction process and plots
    '''
    print('''
    #########################
    #         Test 1        #
    #########################
    ''')
    subnetworks, peak_dict, edge_dict = test_read_url()
    big = [g for g in subnetworks if g.size()>15]

    WV = Weavor(peak_dict, isotope_search_patterns, adduct_search_patterns, 
                mz_tolerance_ppm=5, mode='pos')

    print(WV.mzgrid)

    print("\n\n")
    print("Example khipu of one empirical compound from demo data.")
    KP = khipu_diagnosis(big[ 1 ])
    print("Input network edges: ")
    print(KP.input_network.edges(data=True))

    KP.build_khipu(WV)

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

    print("Feature intensities: ")
    print(KP.get_khipu_intensities())
    # KP.plot_khipu_diagram()

    print("\n\n")
    print("Multiple example khipus: ")
    print("======================== \n")
    for g in big[ -6: ]:
        KP = Khipu(g)
        KP.build_khipu(WV)
        print(KP.sorted_mz_peak_ids, "\n")
        KP.print_khipu()
        print('\n\n')

    print("\n\n")
    print("Build and export JSON file of full list of khipus: ")
    print("================================================== \n")

    khipu_list = graphs_to_khipu_list(
        subnetworks, WV, mz_tolerance_ppm=5)
    
    # khipu_list = extend_khipu_list(khipu_list, peak_dict, adduct_search_patterns_extended)
    print("\n\n ~~~~~~ Got %d khipus ~~~~~~~ \n\n" %len(khipu_list))
    empCpds = export_empCpd_khipu_list(khipu_list)

    outfile = 'khipu_test_empricalCompounds.json'
    with open(outfile, 'w', encoding='utf-8') as f:
        json.dump(empCpds, f, ensure_ascii=False, indent=2)

def test2(mz_tolerance_ppm = 5, mode = 'pos'):
    '''Test epdsConstructor wrapper class
    '''
    print('''
    #########################
    #         Test 2        #
    #########################
    ''')
    adduct_patterns = adduct_search_patterns
    if mode == 'neg':
        adduct_patterns = adduct_search_patterns_neg

    list_peaks = get_peaklist_from_url(url)

    ECCON = epdsConstructor(list_peaks, mode=mode)
    dict_empCpds = ECCON.peaks_to_epdDict(
                        isotope_search_patterns,
                        adduct_patterns,
                        extended_adducts,
                        mz_tolerance_ppm,
                        # rt_tolerance,
        ) 
    outfile = 'khipu_test2.json'
    with open(outfile, 'w', encoding='utf-8') as f:
        json.dump(dict_empCpds, f, ensure_ascii=False, indent=2)


#
# -----------------------------------------------------------------------------
#
if __name__ == '__main__':

    print("\n\n~~~~ testing Khipu (%s) ~~~\n" %__version__)
    test1()
    test2()