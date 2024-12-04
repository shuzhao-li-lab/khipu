# this module contains methods for analyzing khipu data

import numpy as np
import json
import isocor
import pandas as pd
import requests as req

from mass2chem.formula import parse_chemformula_dict, dict_to_hill_formula
from intervaltree import IntervalTree
from .utils import ADDUCT_TO_FORMULA_DELTAS
from .isotopes import get_isotope_data

TRACER_ELEMENT_MAP = {
    "13C": "C",
    "15N": "N",
}

TRACER_ISOTOPOLOGUE_MAP = {
    "13C": "13C/12C",
    "15N": "15N/14N",
}

TRACER_ABUNDANCE_VECTORS = {
    "13C": [0, 100],
    "15N": [0, 100],
}


def detect_labelling(khipu_list, 
                     unlabeled_samples, 
                     labeled_samples, 
                     result_name, 
                     labeling_threshold=10, 
                     skip_isos=None):
    """
    Run detect labelling on every empirical compound.

    Args:
        khipu_list (_type_): _description_
        unlabeled_samples (list): list of unlabeled samples
        labeled_samples (list): list of labeled samples
        result_name (str): save the results using this key
        labeling_threshold (int, optional): the ratio the mean labeled sample intensity must exceed the
        mean unlabeled sample intensity.. Defaults to 10.
        skip_isos (list, optional): isotoplogues to skip in the comparison. Defaults to ["M0"].

    Updates:
        the "labeling_scores" field in every khipu by adding the results under the key "result_name"
    """    
    skip_isos = {"M0"} if skip_isos is None else skip_isos
    for khipu in khipu_list:
        detect_labelling_khipu(khipu, 
                               unlabeled_samples, 
                               labeled_samples, 
                               result_name, 
                               labeling_threshold, 
                               skip_isos)
    return khipu_list

def detect_labelling_khipu(khipu, 
                           unlabeled_samples, 
                           labeled_samples, 
                           result_name, 
                           labeling_threshold=10, 
                           skip_isos=None):
    """
    This counts the number of presumed labeled and unlabeled isotopologue features in each empCpd.

    This compares the mean intensity of unlabeled and labeled samples, if the labeled is greater than
    unlabeled mean intensity * labeling_threhsold, this is a "good" labelled isotopologue, else it
    is a "bad" labelled isotopologue. The results on saved in a dictionary called "labeling_scores"
    using the result_name field with counters for "good_labelling" and "bad_labelling". Isotopologues
    in the skip_isos list or set, will not be considered in the counting. Should always include M0. 

    Args:
        khipu (dict): an empcpd
        unlabeled_samples (list): list of unlabeled samples
        labeled_samples (list): list of labeled samples
        result_name (str): save the results using this key
        labeling_threshold (int, optional): the ratio the mean labeled sample intensity must exceed the
        mean unlabeled sample intensity.. Defaults to 10.
        skip_isos (list, optional): isotoplogues to skip in the comparison. Defaults to ["M0"].

    Updates:
        the "labeling_scores" field in the khipu by adding the results under the key "result_name"
    """    
    skip_isos = {"M0"} if skip_isos is None else skip_isos
    good_isotopologues = 0
    bad_isotopologues = 0
    for peak in khipu["MS1_pseudo_Spectra"]:
        isotope = peak["isotope"]
        if isotope not in skip_isos:
            unlabeled_avg = np.mean([peak[us] for us in unlabeled_samples])
            labeled_avg = np.mean([peak[s] for s in labeled_samples])
            if labeled_avg > unlabeled_avg * labeling_threshold:
                good_isotopologues += 1
            else:
                bad_isotopologues += 1
    if "labeling_scores" not in khipu:
        khipu["labeling_scores"] = {}
    khipu["labeling_scores"][result_name] = {
        "good_labelling": good_isotopologues,
        "bad_labelling": bad_isotopologues
    }

def correct_natural_abundance(khipu_list, 
                              unlabeled_samples, 
                              labeled_samples, 
                              tracer,
                              tracer_purity,
                              resolution,
                              mz_of_resolution,
                              resolution_formula_code,
                              isotope_data_path=None,
                              unique_only=True):
    """
    This performs isocor on every khipu. See the correct_natural_abundance_khipu
    function for more details. 

    Args:
        khipu_list (list): list of empcpds
        unlabeled_samples (list): list of unlabeled sample names
        labeled_samples (list): list of labeled sample names
        tracer (str): the isotope to correct for
        tracer_purity (float): the purity of the tracer (assumes two isotopes per element)
        resolution (integer): resolution of the instrument in resolving power
        mz_of_resolution (float): the m/z at which the resolution was measured
        resolution_formula_code (str): controls the type of instrument can be 'orbitrap', 'ft-icr', or 'constant'
        charge (int): the charge of the molecule
        isotope_data_path (str, optional): if provided, use the isotope data in this file. Defaults to None.
        unique_only (bool): if true, only correct the empcpds with unique formula annotations

    Updates:
        This creates a new set of feature intensities for the labelled samples based on the correction for 
        each adduct, formula. 
    """
    L = []
    for khipu in khipu_list:
        new_khipu = correct_natural_abundance_khipu(khipu, 
                                                    unlabeled_samples, 
                                                    labeled_samples, 
                                                    tracer,
                                                    tracer_purity,
                                                    resolution,
                                                    mz_of_resolution,
                                                    resolution_formula_code,
                                                    isotope_data_path=isotope_data_path,
                                                    unique_only=unique_only,)
        L.append(new_khipu)
    return L

#@functools.lru_cache(1)
def __build_isocor_corrector(formula, 
                             tracer, 
                             tracer_purity, 
                             resolution, 
                             mz_of_resolution,
                             resolution_formula_code,
                             charge,
                             isotope_data_path=None,
                             adduct=None):
    """
    This function returns the appropriate isocor corrector for the given
    input parameters.

    Args:
        formula (str): the formula of the compound to correct
        tracer (str): the isotope to correct for
        tracer_purity (float): the purity of the tracer (assumes two isotopes per element)
        resolution (integer): resolution of the instrument in resolving power
        mz_of_resolution (float): the m/z at which the resolution was measured
        resolution_formula_code (str): controls the type of instrument can be 'orbitrap', 'ft-icr', or 'constant'
        charge (int): the charge of the molecule
        isotope_data_path (str, optional): if provided, use the isotope data in this file. Defaults to None.

    Returns:
        isocor corrector: an instance of an isocor corrector
    """     
    
    return isocor.mscorrectors.MetaboliteCorrectorFactory(
        formula=formula,
        tracer=tracer,
        label=formula,
        inchi=None,
        data_isotopes=get_isotope_data(isotope_data_path),
        derivative_formula=adduct,
        tracer_purity=[1 - tracer_purity, tracer_purity] if tracer_purity else None,
        correct_NA_tracer=True,
        resolution=resolution,
        mz_of_resolution=mz_of_resolution,
        resolution_formula_code=resolution_formula_code,
        charge=charge
    )

def correct_natural_abundance_khipu(khipu, 
                                    unlabeled_samples, 
                                    labeled_samples, 
                                    tracer,
                                    tracer_purity,
                                    resolution,
                                    mz_of_resolution,
                                    resolution_formula_code,
                                    isotope_data_path,
                                    unique_only,):
    
    """
    This performs isocor on the khipu. Only the labeled_samples are corrected.

    The correction is performed using isocor, and the adduct is assumed to be 
    unlabeled. 

    The resolution of the instrument, the mz_of_resolution, and resolution
    formula codes are used for constructing the constructor. If unique only is true
    a khipu is only corrected if the khipu has a single formula annotated to it. 

    Args:
        khipu_list (list): list of empcpds
        unlabeled_samples (list): list of unlabeled sample names
        labeled_samples (list): list of labeled sample names
        tracer (str): the isotope to correct for
        tracer_purity (float): the purity of the tracer (assumes two isotopes per element)
        resolution (integer): resolution of the instrument in resolving power
        mz_of_resolution (float): the m/z at which the resolution was measured
        resolution_formula_code (str): controls the type of instrument can be 'orbitrap', 'ft-icr', or 'constant'
        charge (int): the charge of the molecule
        isotope_data_path (str, optional): if provided, use the isotope data in this file. Defaults to None.
        unique_only (bool): if true, only correct the empcpds with unique formula annotations

    Updates:
        This creates a new set of feature intensities for the labelled samples based on the correction for 
        each adduct, formula. 
    """


    if "isocor_results" not in khipu:
        khipu["isocor_results"] = {}

    if '3x' in khipu["MS1_pseudo_Spectra"][0]:
        charge = 3
    elif '2x' in khipu["MS1_pseudo_Spectra"][0]:
        charge = 2
    else:
        charge = 1

    if "list_matches" in khipu and khipu["list_matches"]: 
        # yuanye, will need to update for your annotation format TODO
        formulas = [x[0].split("_")[0] for x in khipu["list_matches"]]
        if len(formulas) > 1 and unique_only:
            return khipu
        abundance_vectors = {}
        peak_lookup = {x['id']: x for x in khipu["MS1_pseudo_Spectra"]}
        for peak in khipu["MS1_pseudo_Spectra"]:
            adduct = peak["modification"]
            # "M+H+, 3x charged". now it seems in one khipu all features have the same charge. Will change later
            if 'charged' in adduct:
                adduct = adduct.split(',')[0]
            if adduct not in abundance_vectors:
                abundance_vectors[adduct] = {}
            isotope = peak["isotope"]
            if "M0" in isotope:
                count = 0
            elif "*" not in isotope:
                count = 1
            else:
                count = int(isotope.split(",")[0].rstrip().split("*")[-1])
            abundance_vectors[adduct][int(count)] = peak["id"]
        for formula in formulas:
            for adduct, peaks_for_adduct in abundance_vectors.items():
                khipu['isocor_results'][formula + "_" + adduct] = {}
                metabolite_formula_dict = parse_chemformula_dict(formula)
                adduct_formula_dict = {}
                combined_formula_dict = {k: v for k,v in metabolite_formula_dict.items()}
                for ele, count in ADDUCT_TO_FORMULA_DELTAS[adduct][2].items():
                    # if negative, this means it was removed from the molecule being analyzed
                    if count < 0:
                        metabolite_formula_dict[ele] = metabolite_formula_dict.get(ele, 0) + count
                    # if positive, it was added, implying it came from the environment and thus cannot be labeled.
                    else:
                        adduct_formula_dict[ele] = adduct_formula_dict.get(ele, 0) + count
                    combined_formula_dict[ele] = combined_formula_dict.get(ele, 0) + count
                adduct_formula = dict_to_hill_formula(adduct_formula_dict)
                corrected_metabolite_formula = dict_to_hill_formula(metabolite_formula_dict)
                adducted_metabolite_formula = dict_to_hill_formula(combined_formula_dict)
                max_ele = parse_chemformula_dict(adducted_metabolite_formula).get(TRACER_ELEMENT_MAP[tracer], 0)
                if max_ele:
                    corrector = __build_isocor_corrector(adducted_metabolite_formula, 
                                                    tracer,
                                                    tracer_purity,
                                                    resolution,
                                                    mz_of_resolution,
                                                    resolution_formula_code,
                                                    charge,
                                                    isotope_data_path=isotope_data_path,
                                                    adduct=adduct_formula)
                    for ls in labeled_samples:
                        peak_vector = [peaks_for_adduct.get(i, None) for i in range(max_ele + 1)]
                        to_correct = [peak_lookup.get(peaks_for_adduct.get(i, None), {}).get(ls, 0) for i in range(max_ele + 1)]
                        corrected_area, _, _, _ = corrector.correct(to_correct)
                        for i, (corr_intensity, f_id) in enumerate(zip(corrected_area, peak_vector)):
                            if corr_intensity > 0:
                                if i == 0:
                                    iso = "M0"
                                elif i == 1: 
                                    iso = TRACER_ISOTOPOLOGUE_MAP[tracer]
                                else:
                                    iso = TRACER_ISOTOPOLOGUE_MAP[tracer] + "*" + str(i)
                                if charge > 1:
                                    charge_string = ", " + str(charge) + "x charged"
                                    iso += charge_string
                                if iso not in khipu["isocor_results"][formula + "_" + adduct]:
                                    if f_id is not None:
                                        khipu["isocor_results"][formula + "_" + adduct][iso] = dict(peak_lookup[f_id]) # deep copy, otherwise original peak will be changed
                                        khipu["isocor_results"][formula + "_" + adduct][iso]['nac_intensities'] = {}
                                        for sample in labeled_samples:
                                            khipu["isocor_results"][formula + "_" + adduct][iso]['nac_intensities'][sample] = 0
                                        for sample in unlabeled_samples:
                                            if sample in khipu["isocor_results"][formula + "_" + adduct][iso]['nac_intensities'].keys():
                                                del khipu["isocor_results"][formula + "_" + adduct][iso]['nac_intensities'][sample]
                                    else:
                                        khipu["isocor_results"][formula + "_" + adduct][iso] = {
                                            "apex": None,
                                            "peak_area": None,
                                            "height": None,
                                            "left_base": None,
                                            "right_base": None,
                                            "goodness_fitting": None,
                                            "cSelectivity": None,
                                            "mz": None,
                                            "snr": None,
                                            "id_number": "F_" + formula + "_" + iso + "_" + adduct,
                                            "rtime": None,
                                            "rtime_left_base": None,
                                            "rtime_right_base": None,
                                            "id": "F_" + formula + "_" + iso + "_" + adduct,
                                            "isotope": iso,
                                            "modification": adduct,
                                            "ion_relation": iso + "," + adduct,
                                            "parent_epd_id": None
                                        }
                                        khipu["isocor_results"][formula + "_" + adduct][iso]['nac_intensities'] = {}
                                khipu["isocor_results"][formula + "_" + adduct][iso]['nac_intensities'][ls] = corr_intensity
                                khipu["isocor_results"][formula + "_" + adduct][iso][ls] = corr_intensity
                khipu["isocor_results"][formula + "_" + adduct] = list(khipu["isocor_results"][formula + "_" + adduct].values())
    return khipu
     
     
def measure_overlap_with_GEM(khipu_list, GEM, mz_tol=5):
    """
    This function performs annotation against the provided GEM
    using IntervalTrees with the inferred neutral mass and 
    the specified mz_tolerance.

    Args:
        khipu_list (list): list of empirical compounds
        GEM (str): path to GEM in JSON format
        mz_tol (int, optional): mz tolerance for the query. Defaults to 5.

    Updates:
        "in_GEM" and "GEM_annotations" a boolean and list, in the 
        empCpd.
    """    
    GEM_data = json.load(open(GEM))
    GEM_mz_tree = IntervalTree()
    id_to_name = {}
    for cpd in GEM_data["list_of_compounds"]:
        id_to_name[cpd['id']] = cpd['name']
        mass = cpd["neutral_mono_mass"]
        if mass:
            mass_err = abs(mass / 1e6 * mz_tol)
            GEM_mz_tree.addi(mass - mass_err, mass + mass_err, cpd['id'])
    for khipu in khipu_list:
        if "in_GEM" not in khipu:
            khipu["in_GEM"] = False
        for r in GEM_mz_tree.at(khipu["neutral_formula_mass"]):
            if "GEM_annotations" not in khipu:
                khipu["GEM_annotations"] = []
            khipu["in_GEM"] = True
            khipu["GEM_annotations"].append(id_to_name[r.data])
    return khipu_list
