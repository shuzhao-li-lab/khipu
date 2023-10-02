'''
Ad hoc function to clean up current khipu or asari output, 
1. add neutral formula to empCpd dictionary
2. check possible isotope pattern against number of atoms in the formula and add tag 'okay_formula'

Will add mass calibration step too.

These will be further developed and added to future khipu versions.

SL 2023-10-01
'''

from mass2chem.search import find_best_match_centurion_indexed_list

from khipu.epdsConstructor import epdsConstructor
from khipu.utils import adduct_search_patterns, isotope_search_patterns, extended_adducts
from khipu.plot import plot_json_khipu
from khipu.tree_formulas import tree_hmdb4, tree_pubchemlite


def get_khipu_dict(features, mode='pos'):
    ECON = epdsConstructor(features, mode=mode)
    khipu_dict = ECON.peaks_to_epdDict(
        isotope_search_patterns = isotope_search_patterns,
        adduct_search_patterns = adduct_search_patterns,
        extended_adducts = extended_adducts,
        mz_tolerance_ppm=5, 
        rt_tolerance=1, 
        charges=[1, 2, 3],
        has_parent_masstrack=False,
    )
    return khipu_dict


def get_M0(MS1_pseudo_Spectra):
    '''returns M0 feature with highest representative_intensity.
    Without verifying which ion form.'''
    M0 = [(f['representative_intensity'], f['id'], f) for f in 
          MS1_pseudo_Spectra if f['isotope']=='M0']
    if M0:    # use f['id'] as tie breaker in sort
        return sorted(M0, reverse=True)[0][2]
    else:
        return []
    
def get_M1(MS1_pseudo_Spectra):
    '''returns M+1 feature with highest representative_intensity.
    Without verifying which ion form.'''
    M = [(f['representative_intensity'], f['id'], f) for f in 
          MS1_pseudo_Spectra if f['isotope']=='13C/12C']
    if M:
        return sorted(M, reverse=True)[0][2]
    else:
        return []
     
def get_highest_13C(MS1_pseudo_Spectra):
    '''returns 13C labeled feature with highest representative_intensity.
    Without verifying which ion form. Because the label goes with sepecific atoms depending on pathway.
    '''
    M = [(f['representative_intensity'], f) for f in 
          MS1_pseudo_Spectra if '13C/12C' in f['isotope']]
    if M:
        return sorted(M, reverse=True)[0][1]
    else:
        return []

def get_labeled(kdict, label_ratio_filter=0.2):
    '''
    kdict : khipu_dict as input
    returns 
    list of khipus with 13C/12C ratio greater than label_ratio_filter, using representative_intensity.
    '''
    khipus_good_ratios = []
    for k,v in kdict.items():
        # interim_id = v['interim_id']
        M0, Mx = get_M0(v['MS1_pseudo_Spectra']), get_highest_13C(v['MS1_pseudo_Spectra'])
        if M0 and Mx:
            try:
                ratio = Mx['representative_intensity']/(M0['representative_intensity'] + 0.1)
                if ratio > label_ratio_filter:
                    khipus_good_ratios.append( (k, ratio) )
            except IndexError:
                pass
    print(len(khipus_good_ratios))
    return khipus_good_ratios
    
def filter_khipus(kdict, unlabelled=[1, 2, 3], labeled=[0, 4, 5], natural_ratio_limit=0.5):
    '''
    kdict : khipu_dict as input
    
    returns 
    list of khipus with good natural 13C ratio, list of khipus with increased 13C ratio.
    '''
    khipus_good_natural_ratios, khipus_increased_ratios = [], []
    for k,v in kdict.items():
        # interim_id = v['interim_id']
        M0, M1, Mx = get_M0(v['MS1_pseudo_Spectra']), get_M1(v['MS1_pseudo_Spectra']
                                            ), get_highest_13C(v['MS1_pseudo_Spectra'])
        if M0 and M1:
            try:
                unlabelled_13C = np.array(M1['intensities'])[unlabelled].mean()
                unlabelled_12C = np.array(M0['intensities'])[unlabelled].mean()
                ratio_natural = unlabelled_13C/(unlabelled_12C+0.1)
                if ratio_natural < natural_ratio_limit:
                    khipus_good_natural_ratios.append( (k, ratio_natural) )
                    if Mx:
                        labelled_13C = np.array(Mx['intensities'])[labeled].mean()
                        labelled_12C = np.array(M0['intensities'])[labeled].mean()
                        ratio_labeled = labelled_13C/(labelled_12C+0.1)
                        if ratio_labeled > ratio_natural:
                            khipus_increased_ratios.append( (k, ratio_labeled) )
            except IndexError:
                pass
    print(len(khipus_good_natural_ratios), len(khipus_increased_ratios))
    return khipus_good_natural_ratios, khipus_increased_ratios


def find_formula(khipu_dict, tree_hmdb4, tree_pubchemlite, limit_ppm=10):
    '''Find a matched formula from tree_hmdb4 or tree_pubchemlite
    '''
    for k, v in khipu_dict.items():
        match_hmdb4 = find_best_match_centurion_indexed_list(v['neutral_formula_mass'], tree_hmdb4, limit_ppm)
        if match_hmdb4:
            # expecting e.g. {'mz': 222.073953, 'id_number': 'C8H14O7_222.073953'}
            v['neutral_formula'] = match_hmdb4['id_number'].split('_')[0]
            v['Database_referred'].append('hmdb4')
            v['identity'].append(match_hmdb4['id_number'])
        else:
            match_pcl = find_best_match_centurion_indexed_list(v['neutral_formula_mass'], tree_pubchemlite, limit_ppm)
            if match_pcl:
                v['neutral_formula'] = match_pcl['id_number'].split('_')[0]
                v['Database_referred'].append('PubChemLite')
                v['identity'].append(match_pcl['id_number'])
                
    return khipu_dict

