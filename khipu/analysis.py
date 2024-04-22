# this module contains methods for analyzing khipu data

import numpy as np
import json
import isocor
import pandas as pd
import requests as req
import re
import os
from decimal import Decimal
from mass2chem.formula import parse_chemformula_dict, dict_to_hill_formula
from intervaltree import IntervalTree
from .utils import ADDUCT_TO_FORMULA_DELTAS

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

def get_all_isotope_data(nist_path='.'):
    def _stripColNames(df):
        df.rename(columns=lambda x: x.strip())

    def _stripCol(df, listcolumns):
        for col in listcolumns:
            df[col] = df[col].str.strip()

    def _makeIsotopesDict(df):
        gb = df.groupby('element')
        d = {}
        for g in gb.groups:
            dg = gb.get_group(g).sort_values(by='mass').to_dict('list')
            e = dg.pop('element')[0]
            d[e] = dg
        return d
    
    with open(str(nist_path), 'r', encoding='utf-8') as fp:
        dfIsotopes = pd.read_csv(fp, converters={'mass': Decimal, 'abundance': np.float64})
    _stripColNames(dfIsotopes)
    _stripCol(dfIsotopes, ['element', ])
    dictIsotopes = _makeIsotopesDict(dfIsotopes)
    return dictIsotopes

def detect_labelling(khipu_list, unlabeled_samples, labeled_samples, result_name, labeling_threshold=10, skip_isos=None):
    skip_isos = {"M0"} if skip_isos is None else skip_isos
    for khipu in khipu_list:
        detect_labelling_khipu(khipu, unlabeled_samples, labeled_samples, result_name, labeling_threshold, skip_isos)
    return khipu_list

def detect_labelling_khipu(khipu, unlabeled_samples, labeled_samples, result_name, labeling_threshold=10, skip_isos=["M0"]):
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
                              nist_path,
                              unique_only):
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
                                                    nist_path,
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
                             nist_path):
     return isocor.mscorrectors.MetaboliteCorrectorFactory(
        formula=formula,
        tracer=tracer,
        label=formula,
        inchi=None,
        data_isotopes=get_all_isotope_data(nist_path),
        derivative_formula=None,
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
                                    nist_path,
                                    unique_only,):
    if "isocor_results" not in khipu:
        khipu["isocor_results"] = {}

    if '3x' in khipu["MS1_pseudo_Spectra"][0]:
        charge = 3
    elif '2x' in khipu["MS1_pseudo_Spectra"][0]:
        charge = 2
    else:
        charge = 1

    if "list_matches" in khipu and khipu["list_matches"]: # yuanye, will need to update
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
                formula_dict = parse_chemformula_dict(formula)

                # the formula should have the element of the tracer
                tracer_ele = re.findall(r'[A-Z]+', tracer)
                if not tracer_ele[0] in formula_dict.keys():
                    continue

                for ele, count in ADDUCT_TO_FORMULA_DELTAS[adduct][2].items():
                    formula_dict[ele] = formula_dict.get(ele, 0) + count
                adduct_corrected_formula = dict_to_hill_formula(formula_dict)
                corrector = __build_isocor_corrector(adduct_corrected_formula, 
                                                    tracer,
                                                    tracer_purity,
                                                    resolution,
                                                    mz_of_resolution,
                                                    resolution_formula_code,
                                                    charge,
                                                    nist_path)
                max_ele = parse_chemformula_dict(adduct_corrected_formula).get(TRACER_ELEMENT_MAP[tracer], 0)
                if max_ele:
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
                khipu["isocor_results"][formula + "_" + adduct] = list(khipu["isocor_results"][formula + "_" + adduct].values())
    return khipu
     
     
def measure_overlap_with_GSMM(khipu_list, GEM, mz_tol=5):
    GEM_data = json.load(open(GEM))
    GEM_mz_tree = IntervalTree()
    for cpd in GEM_data["list_of_compounds"]:
        mass = cpd["neutral_mono_mass"]
        if mass:
            mass_err = mass / 1e6 * mz_tol
            GEM_mz_tree.addi(mass - mass_err, mass + mass_err, cpd['id'])
    for khipu in khipu_list:
        khipu["in_GEM"] = bool(GEM_mz_tree.at(khipu["neutral_formula_mass"]))
    return khipu_list
