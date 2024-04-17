# this module contains methods for analyzing khipu data

import numpy as np
import json
import isocor
import functools
from mass2chem.formula import parse_chemformula_dict
from intervaltree import IntervalTree

def detect_labelling(khipu_list, unlabeled_samples, labeled_samples, result_name, labeling_threshold=10, skip_isos=["M0"]):
    L = []
    for khipu in khipu_list:
        new_khipu = detect_labelling_khipu(khipu, unlabeled_samples, labeled_samples, result_name, labeling_threshold, skip_isos)
        L.append(new_khipu)
    return L

def detect_labelling_khipu(khipu, unlabeled_samples, labeled_samples, result_name, labeling_threshold=10, skip_isos=["M0"]):
    skip_iso_set = set(skip_isos)
    good_isotopologues = 0
    bad_isotopologues = 0
    for peak in khipu["MS1_pseudo_Spectra"]:
        ion_relation = peak["ion_relation"]
        isotope, _ = ion_relation.split(",")
        if isotope not in skip_iso_set:
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
    return khipu

def correct_natural_abundance(khipu_list, 
                              unlabeled_samples, 
                              labeled_samples, 
                              tracer,
                              tracer_purity,
                              resolution,
                              mz_of_resolution,
                              resolution_formula_code,
                              ):
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
                                                    unique_only=False)
        L.append(new_khipu)
    return L

#@functools.lru_cache(1)
def __build_isocor_corrector(formula, 
                             tracer, 
                             tracer_purity, 
                             resolution, 
                             mz_of_resolution,
                             resolution_formula_code,
                             charge):
    return isocor.mscorrectors.MetaboliteCorrectorFactory(
        formula=formula,
        tracer=tracer,
        label=formula,
        inchi=None,
        data_isotopes=None,
        derivative_formula=None,
        tracer_purity=[1 - tracer_purity, tracer_purity] if tracer_purity else None,
        correct_NA_tracer=False,
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
                                    unique_only=False):
    if '3x' in khipu["MS1_pseudo_Spectra"][0]:
        charge = 3
    elif '2x' in khipu["MS1_pseudo_Spectra"][0]:
        charge = 2
    else:
        charge = 1

    if "list_matches" in khipu and khipu["list_matches"]:
        formulas = [x[0].split("_")[0] for x in khipu["list_matches"]]
        if len(formulas) > 1 and unique_only:
            return khipu
        abundance_vectors = {}
        peak_lookup = {x['id']: x for x in khipu["MS1_pseudo_Spectra"]}
        for peak in khipu["MS1_pseudo_Spectra"]:
            adduct = peak["modification"]
            if adduct not in abundance_vectors:
                abundance_vectors[adduct] = {}
            isotope = peak["isotope"]
            if "M0" in isotope:
                count = 0
            elif "13C/12C" in isotope:
                count = 1
            else:
                count = int(isotope.split(",").rstrip().split("*")[-1])
            abundance_vectors[adduct][int(count)] = peak["id"]
        for formula in formulas:
            if "corrected_for_" + formula not in khipu:
                khipu["corrected_for_" + formula] = {}
            corrector = __build_isocor_corrector(formula, 
                                                tracer,
                                                tracer_purity,
                                                resolution,
                                                mz_of_resolution,
                                                resolution_formula_code,
                                                charge)
            max_C = parse_chemformula_dict(formula).get("C", 0)
        
            if max_C:
                for ls in labeled_samples:
                    for adduct, peaks_for_adduct in abundance_vectors.items():
                        correction_vector = []
                        peak_vector = []
                        for i in range(max_C + 1):
                            f_id = peaks_for_adduct.get(i, None)
                            peak_vector.append(f_id)
                            correction_vector.append(peak_lookup.get(f_id, {}).get(ls, 0))
                        corrected_area, _, _, _ = corrector.correct(correction_vector)
                        if adduct not in khipu["corrected_for_" + formula]:
                            khipu["corrected_for_" + formula][adduct] = {}
                        #print(correction_vector)
                        #print(peak_vector)
                        for i, (corr_intensity, f_id) in enumerate(zip(corrected_area, peak_vector)):
                            #print(i, corr_intensity, f_id)
                            if corr_intensity > 0:
                                if i == 0:
                                    iso = "M0"
                                elif i == 1: 
                                    iso = "13C/12C"
                                else:
                                    iso = "13C/12C*" + str(i)
                                if charge > 1:
                                    charge_string = ", " + str(charge) + "x charged"
                                    iso += charge_string
                                if iso not in khipu["corrected_for_" + formula][adduct]:
                                    khipu["corrected_for_" + formula][adduct][iso] = {}
                                if "original_f_id" not in khipu["corrected_for_" + formula][adduct][iso]:
                                    khipu["corrected_for_" + formula][adduct][iso]["original_f_id"] = f_id
                                    khipu["corrected_for_" + formula][adduct][iso]["adduct"] = adduct
                                    khipu["corrected_for_" + formula][adduct][iso]["isotope"] = iso
                                khipu["corrected_for_" + formula][adduct][iso][ls] = corr_intensity
                                #print("***", adduct, iso, corr_intensity)
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
        khipu["in_GEM"] = bool(GEM_mz_tree.at("neutral_formula_mass"))
    return khipu_list


