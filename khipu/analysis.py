# this module contains methods for analyzing khipu data

import numpy as np
import json
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
            