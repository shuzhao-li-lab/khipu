from intervaltree import IntervalTree
from itertools import product
import numpy as np
import pandas as pd
from copy import deepcopy
from mass2chem.formula import parse_chemformula_dict
import json
import uuid
import tqdm
import sys

class MassDeltaExplainer(object):
    def __init__(self, 
                 component_path, 
                 PPM_MZ_TOL=5, 
                 MIN_ISO_PROB=0.001, 
                 ELEMENT_NEUTRAL=False, 
                 MASS_LIMIT=200, 
                 DEPTH=3, 
                 CLEAN_SOLUTIONS=True,
                 FIXED_MZ_ERR=0.0001) -> None:
        self.components = pd.read_csv(component_path).to_dict(orient='records')
        for comp in self.components:
            comp['id'] = str(uuid.uuid4())
        self.mz_tree = IntervalTree() # this is used for querying mass deltas efficiently
        self.mz_tree_depth = DEPTH # the max number of components that may be combined to form a solution
        self.PPM_MZ_TOL = PPM_MZ_TOL # the mz tolerance used for searching, in ppm
        self.MIN_ISO_PROB = MIN_ISO_PROB # solutions whose probability product (this is not a precise measurement) is below this value are discarded
        self.ELEMENT_NEUTRAL = ELEMENT_NEUTRAL # if true, consider only solutions in which the net formula did not change
        self.MASS_LIMIT = MASS_LIMIT # eliminate solutions whose mass delta exceeds this value
        self.CLEAN_SOLUTIONS = CLEAN_SOLUTIONS # remove bookkeeping fields from output
        self.FIXED_MZ_ERR = FIXED_MZ_ERR
        self.combine_components() # build the mz_tree

    def combine_components(self):
        self.mz_tree = IntervalTree()
        
        # visualize a tree with depth+1 layers of nodes, let the root node be empty.
        # permuting all combinations of components is thus a breadth first search of the tree starting at the root.
        # each path from root to leaf is a solution, or combination of components to be considered.
        # the code below performs this search and returns all valid solutions based on the rules set in __init__

        working_depth = 0
        all_solutions = [{
            "mass": 0, # this mass is for generation and will be recomputed to prevent floating point issues
            "components": {
                "added": [],
                "removed": []
            },
            "working_prob": 1.0 # this prob is for generation, will be recomputed to prevent floating point issues
        }]
        while working_depth < self.mz_tree_depth:
            for solution in tqdm.tqdm(list(all_solutions), desc="enumerating depth " + str(working_depth)):
                for component in self.components:
                    for mode in ["added", "removed"]:
                        new_solution = deepcopy(solution)
                        new_solution["mass"] += component["mz_delta"] * (1 if mode == "added" else -1)
                        new_solution["working_prob"] *= component["prob"]
                        new_solution["components"][mode].append(component)
                        # filter any degenerate solution, in which, we have added and removed the same component
                        if not set(x['id'] for x in new_solution["components"]['added']).intersection(set(x['id'] for x in new_solution["components"]['removed'])):
                            if new_solution["working_prob"] > self.MIN_ISO_PROB and new_solution["mass"] < self.MASS_LIMIT:
                                all_solutions.append(new_solution)
            working_depth += 1
        
        added = set()
        working_solutions = []
        for solution in all_solutions:
            #print(json.dumps(solution, indent=4))
            solution["solution_mass_delta"] = np.sum([x['mz_delta'] for x in solution["components"]["added"]] + [-1 * x['mz_delta'] for x in solution["components"]["removed"]])
            solution["solution_prob"] = np.prod([x['prob'] for x in solution["components"]["added"]] + [x['prob'] for x in solution["components"]["removed"]])
            solution["isotope_only"] = self.isotope_only_solution(solution)
            #print(json.dumps(solution, indent=4))

            if self.CLEAN_SOLUTIONS:
                del solution['working_prob']
                del solution['mass']

            if solution["solution_prob"] > self.MIN_ISO_PROB: # filter improbable
                if (self.ELEMENT_NEUTRAL and solution['isotope_only']) or not self.ELEMENT_NEUTRAL: # filter modifications that change formulas
                    mz_err = abs(solution['solution_mass_delta'] / 1e6 * self.PPM_MZ_TOL) # make interval based on ppm
                    mz_err += self.FIXED_MZ_ERR
                    if mz_err > 0: # to avoid null intervals
                        comps_added = sorted([x["Name"] for x in solution['components']["added"]])
                        comps_removed = sorted([x["Name"] for x in solution['components']["removed"]])
                        key = ''
                        if comps_added:
                            key += "added_" + "_".join(comps_added)
                        if comps_added and comps_removed:
                            key += "_"
                        if comps_removed:
                            key += "removed_" + "_".join(comps_removed)
                        if key not in added:
                            added.add(key)
                            solution["key"] = key
                            self.mz_tree.addi(solution['solution_mass_delta'] - mz_err, solution['solution_mass_delta'] + mz_err, solution)
                            working_solutions.append(solution)
        json.dump(working_solutions, open("solutions.json", 'w'), indent=4)

    @staticmethod
    def isotope_only_solution(solution):
        # this checks that the net modification left the elemental formula unchanged
        net_formula = {"added": {}, "removed": {}}
        for mode in ["added", "removed"]:
            for x in solution["components"][mode]:
                for k, v in parse_chemformula_dict(x['formula']).items():
                    if k not in net_formula[mode]:
                        net_formula[mode][k] = 0
                    net_formula[mode][k] += v
        iso_only = True
        for mode in ["added", "removed"]:
            for k in net_formula[mode].keys():
                if net_formula["removed"].get(k, 0) != net_formula["added"].get(k, 0):
                    # at least one element count has differed, thus, not formula neutral
                    iso_only = False
                    break
        return iso_only

    def explains(self, query_mass):
        # given a query mass, find all matching solutions in the tree
        matches = self.mz_tree.at(query_mass)
        return {
            "explained": bool(matches),
            "num_solutions": len(matches),
            "solutions": [x.data for x in matches],
        }
if __name__ == '__main__':
    # example usage
    MDE = MassDeltaExplainer("./components_pos.csv")

    new_ft = []
    for x in pd.read_csv(sys.argv[1], sep="\t").to_dict(orient='records'):
        r1 = MDE.explains(x['delta_mz'])
        #print("\n")
        #print("----")
        #print("SEARCH: ", x['delta_mz'])
        #print("POS:")
        #for z in r1['solutions']:
        #    print("\t", z)

        r2 = MDE.explains(-1 * x['delta_mz'])
        #print("NEG:")
        #for z in r2["solutions"]:
        #    print("\t", z)
        #print(len(r1), len(r2))
        r_combined = {
            "explained": r1['explained'] or r2['explained'],
            "num_solutions": len(r1['solutions']) + len(r2['solutions']),
            "solutions": r1["solutions"] + r2["solutions"]
        }
        if r_combined['explained']:
            x['explained'] = []
            for solution in sorted(r_combined['solutions'], key=lambda x: -x["solution_prob"]):
                net_formula_delta = {}
                for component in solution["components"]["added"]:
                    if component['Name'] not in net_formula_delta:
                        net_formula_delta[component['Name']] = 0
                    net_formula_delta[component['Name']] += 1
                for component in solution["components"]["removed"]:
                    if component['Name'] not in net_formula_delta:
                        net_formula_delta[component['Name']] = 0
                    net_formula_delta[component['Name']] -= 1
                x['explained'].append((float(solution['solution_mass_delta']), float(solution['solution_prob']), net_formula_delta))
            #x['explained'] = [x['key'] for x in sorted(r_combined['solutions'], key=lambda x: -x["solution_prob"])]
        else:
            x['explained'] = []
        new_ft.append(x)

    count_explained = 0
    count_unexplained = 0
    for z in new_ft:
        if z['explained']:
            count_explained += int(z['count_estimate'])
        else:
            count_unexplained += int(z['count_estimate'])


    print(count_explained / (count_explained + count_unexplained))
    print(pd.DataFrame(new_ft).to_csv(sys.argv[1].replace(".tsv", "_results.tsv"), sep="\t"))