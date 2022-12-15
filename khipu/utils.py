'''
Utility functions and m/z patterns.
Adduct rules can be still better learned in the future.
'''
import json
import numpy as np
import networkx as nx
from mass2chem.search import build_centurion_tree, find_all_matches_centurion_indexed_list

adduct_search_patterns = [  (21.9820, 'Na/H'), 
                            (41.026549, 'Acetonitrile'),
                            (17.02655, 'NH3'),
                            (35.9767, 'HCl'),
                            (37.955882, 'K/H'),
                            ]

adduct_search_patterns_neg = [ (34.969402, 'Cl-'), 
                            (44.998201, 'COOH-'),
                            (17.02655, 'NH3'),
                            (21.9820, 'Na/H'), 
                            (41.026549, 'Acetonitrile'),
                            (35.9767, 'HCl'),
                            (37.955882, 'K/H'),
                            ]

isotope_search_patterns = [ (1.003355, '13C/12C', (0, 0.8)),
                            (2.00671, '13C/12C*2', (0, 0.8)),
                            (3.010065, '13C/12C*3', (0, 0.8)),
                            (4.01342, '13C/12C*4', (0, 0.8)),
                            (5.016775, '13C/12C*5', (0, 0.8)),
                            (6.02013, '13C/12C*6', (0, 0.8)),
                            (7.023485, '13C/12C*7', (0, 0.8)),
                            (8.02684, '13C/12C*8', (0, 0.8)),
                            (9.030195, '13C/12C*9', (0, 0.8)),
                            (10.03355, '13C/12C*10', (0, 0.8)),
                            (11.036905, '13C/12C*11', (0, 0.8)),
                            (12.04026, '13C/12C*12', (0, 0.8)),
                            ]

extended_adducts = [(1.0078, 'H'),
                            (10.991, 'Na/H, double charged'),
                            (0.5017, '13C/12C, double charged'),
                            (18.0106, 'H2O'),      # easy to confuse with bio reactions
                            (18.033823, 'M+NH4'),
                            (27.01089904, 'HCN'),
                            (37.94694, 'Ca/H2'),
                            (32.026215, 'MeOH'),
                            (43.96389, 'Na2/H2'),
                            (46.00548, 'HCOOH'),
                            (67.987424, 'NaCOOH'),
                            (83.961361, 'KCOOH'),
                            (97.96737927, 'H2SO4'),
                            (97.97689507, 'H3PO4'),
]


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

def make_peak_tag(peak):
    '''
    peak format: {'id': 'F1', 'mz': 60.0808, 'rtime': 117.7, ...,
    'intensities': [250346.0], 'representative_intensity': 250346.0}
    '''
    return str(round(peak['mz'], 4)) + '@' + str(round(peak['rtime'], 1))

def make_edge_tag(edge):
    '''
    Concatenate sorted str edges by underscore
    '''
    if edge[0] > edge[1]:
        return edge[1] + '_' + edge[0]
    else:
        return edge[0] + '_' + edge[1]

def make_peak_dict(peak_list):
    '''
    Same as search.build_peak_id_dict but uses either 'id' or 'id_number'.
    '''
    if 'id' in peak_list[0]:
        k = 'id'
    elif 'id_number' in peak_list[0]:
        k = 'id_number'
    else:
        raise KeyError("Peaks need 'id' or 'id_number' as key.")
        
    peak_dict = {}
    for p in peak_list:
        p['id'] = p[k]
        peak_dict[p[k]] = p
    return peak_dict

def rt_matched_by_tolerance(P1, P2, rt_tolerance):
    return abs(P2['rtime']-P1['rtime']) < rt_tolerance

def rt_compared_by_values(P1, P2, rt_tolerance=None):
    return P2['rtime'] > P1['rtime']

def assign_masstrack_ids_in_khipu(feature_dict, mz_tolerance_ppm=5):
    '''Assign mass track ids if they are not included in peak dict.
    They should be if features were processed by asari.

    Parameters
    ----------
    feature_dict : feature dictionary indexed by feature ids, based on input network to khipu, 
        thus very limited size.
    mz_tolerance_ppm : ppm tolerance in examining m/z groups.

    Returns
    -------
    Updated feature_dict with 'parent_masstrack_id'.

    Notes
    -----
    Sort by m/z; separate by mz_tolerance_ppm.
    m/z has to be check by ppm because 1) minor variation may exist btw peaks; and 
    2) float numbers are bad for dictionary keys.
    This method can be occasionally problematic in data that were not processed well,
    by accidently merging m/z regions, therefore causing problems in downstream determination of ion relations.
    (Why mass track is good in asari.)
    '''
    features = sorted([(feature_dict[n]['mz'], n) for n in feature_dict])
    _N = len(features)
    current_mz = features[0][0]
    str_mz = str(current_mz)
    feature_dict[features[0][1]]['parent_masstrack_id'] = str_mz
    for x in features[1:]:
        if x[0]-current_mz > current_mz * mz_tolerance_ppm * 0.000001:
            current_mz = x[0]
            str_mz = str(current_mz)
        feature_dict[x[1]]['parent_masstrack_id'] = str_mz

    return feature_dict


def get_isotopic_edge_pairs(list_peaks, 
                    mztree, 
                    search_patterns = [(1.003355, '13C/12C', (0, 0.8))],
                    mz_tolerance_ppm=5, 
                    isotope_rt_tolerance=2, 
                    check_isotope_ratio = False,
                    ):
    '''
    To find all isotope pairs. 
    Similar to search.get_seed_empCpd_signatures, but return unidirectional pairs only.
    If input peaks have overlaps/duplicates, result will contain the redundant overlap peaks.

    Input
    =====
    list_peaks: [{'parent_masstrace_id': 1670, 'mz': 133.09702315984987, 'rtime': 654, 'height': 14388.0, 'id': 555}, ...]
    mztree: indexed list_peaks
    mz_tolerance_ppm: ppm tolerance in examining m/z patterns.
    search_patterns: a list in the format of [(mz difference, notion, (ratio low limit, ratio high limit)), ..]
            This can be obtained through search.isotopic_patterns. The ratios are optional, because 
            1) naturally occuring constrains are based on chemical formula;
            2) rules are different when isotope tracers are introduced to the experiments.
            But it's important to have a comprehensive list here for isotope tracing experiments.
    rt_tolerance: tolerance threshold for deviation in retetion time, arbitrary unit depending on input data.
            Default intended as 2 seconds.

    Return
    ======
    list of lists of peak pairs that match search_patterns patterns, e.g.
    [ (195, 206, '13C/12C'), ...]. 
    '''
    signatures = []
    for P1 in list_peaks:
        matched = [  ] 
        for _pair in search_patterns:
            (mass_difference, relation) = _pair[:2]
            tmp = find_all_matches_centurion_indexed_list(P1['mz'] + mass_difference, mztree, mz_tolerance_ppm)
            for P2 in tmp:
                if abs(P1['rtime']-P2['rtime']) <= isotope_rt_tolerance:
                    if check_isotope_ratio and len(_pair) > 2:  # checking abundance ratio
                        (abundance_ratio_min, abundance_ratio_max) = _pair[2]
                        if abundance_ratio_min*P1['height'] < P2['height'] < abundance_ratio_max*P1['height']:
                            matched.append( (P1['id'], P2['id'], relation) )
                    else:
                        matched.append( (P1['id'], P2['id'], relation) )
        signatures += matched
    return signatures

def get_adduct_edge_pairs(list_peaks, 
                    mztree, 
                    search_patterns = [(1.0078, 'H'), (21.9820, 'Na/H'), (41.026549, 'Acetonitrile')],
                    mz_tolerance_ppm=5, rt_tolerance=2,
                    ):
    '''
    To find all pairs of adducts (fragments and neutral loss can be accommodated using negative mz difference). 

    Input
    =====
    list_peaks: [{'parent_masstrace_id': 1670, 'mz': 133.09702315984987, 'rtime': 654, 'height': 14388.0, 'id': 555}, ...]
    mztree: indexed list_peaks
    mz_tolerance_ppm: ppm tolerance in examining m/z patterns.
    search_patterns: a list in the format of [(mz difference, notion), ...] 
    rt_tolerance: tolerance threshold for deviation in retetion time, arbitrary unit depending on input data.
            Default intended as 2 seconds.

    Return
    ======
    list of lists of peak pairs that match search_patterns patterns, e.g.
    [ (195, 206, 'H/Na'), ...]. 
    '''
    signatures = []
    for P1 in list_peaks:
        matched = [  ] 
        for _pair in search_patterns:
            (mass_difference, relation) = _pair[:2]
            tmp = find_all_matches_centurion_indexed_list(P1['mz'] + mass_difference, mztree, mz_tolerance_ppm)
            for P2 in tmp:
                if abs(P1['rtime']-P2['rtime']) <= rt_tolerance:
                    matched.append( (P1['id'], P2['id'], relation) )
        signatures += matched
    return signatures


def peaks_to_networks(peak_list, 
                    isotope_search_patterns = [ (1.003355, '13C/12C', (0, 0.8)),
                                        (2.00671, '13C/12C*2', (0, 0.8)),
                                        (3.010065, '13C/12C*3', (0, 0.8)),
                                        (4.01342, '13C/12C*4', (0, 0.8)),
                                        (5.016775, '13C/12C*5', (0, 0.8)),
                                        (6.02013, '13C/12C*6', (0, 0.8)),], 
                    adduct_search_patterns = [ (1.0078, 'H'), 
                                        (21.9820, 'Na/H'), 
                                        (41.026549, 'Acetonitrile')
                                        ],
                    mz_tolerance_ppm=5, 
                    rt_tolerance=2, 
                    ):
    '''Search peak_list for patterns of isotopes and adducts, form a network and get connected subnetworks.
    
    Parameters
    ----------
    list_peaks : [{'mz': 133.09702315984987, 'rtime': 654, 'id': 555}, ...]

    isotope_search_patterns : exact list used to retrieve the subnetworks. E.g. 
        [ (1.003355, '13C/12C', (0, 0.8)),
        (2.00671, '13C/12C*2', (0, 0.8)),
        (3.010065, '13C/12C*3', (0, 0.8)),
        (4.01342, '13C/12C*4', (0, 0.8)),
        (5.016775, '13C/12C*5', (0, 0.8)),
        (6.02013, '13C/12C*6', (0, 0.8)),]

    adduct_search_patterns : exact list used to retrieve the subnetworks. 
        It's not recommended to have a long list here, as it's better to search additional 
        in-source modifications after empCpds are seeded. Example adduct_search_patterns list: 
        [ (1.0078, 'H'), 
        (21.9820, 'Na/H'), 
        (41.026549, 'Acetonitrile')]

    mz_tolerance_ppm : ppm tolerance in examining m/z patterns.
    rt_tolerance: tolerance threshold for deviation in retetion time, arbitrary unit depending on input data.
            Default intended as 2 seconds.

    Returns
    -------
    subnetwork : undirected graph. Example edges:
        [('F1606', 'F20', {'type': 'modification', 'tag': 'H'}), 
        ('F3533', 'F20', {'type': 'modification', 'tag': 'Na/H'}), 
        ('F195', 'F20', {'type': 'modification', 'tag': 'Acetonitrile'}), 
        ('F20', 'F807', {'type': 'modification', 'tag': 'Acetonitrile'}), 
        ('F20', 'F53', {'type': 'isotope', 'tag': '13C/12C'}), 
        ('F874', 'F808', {'type': 'isotope', 'tag': '13C/12C'})]

    peak_dict : JSON peaks indexed by ID
    edge_dict : edge_tag is str sorted, but the dict values preserve the direction, which is missed in nx.subnetwork.

    Notes
    -----
    Features of low abundance may not have detectable isotopes, but can have multiple adducts.
    Do not include too many adducts in the intial search. 
    Do not include neutral loss and fragments in initial search.
    They are better done after a list of khipus are constructed.
    '''
    mztree = build_centurion_tree(peak_list)
    peak_dict = make_peak_dict(peak_list)
    iso_edges = get_isotopic_edge_pairs(peak_list, mztree, 
                                        search_patterns=isotope_search_patterns,
                                        mz_tolerance_ppm=mz_tolerance_ppm,
                                        isotope_rt_tolerance=rt_tolerance,
                                        check_isotope_ratio=False,
                                        )
    adduct_edges = get_adduct_edge_pairs(peak_list, mztree,
                                        search_patterns=adduct_search_patterns,
                                        mz_tolerance_ppm=mz_tolerance_ppm,
                                        rt_tolerance=rt_tolerance,
                                        )
    edge_dict = {}
    for e in iso_edges:
        edge_dict[make_edge_tag(e)] = {"source": e[0], "target": e[1], "type": "isotope", "tag": e[2]}
    for e in adduct_edges:
        edge_dict[make_edge_tag(e)] = {"source": e[0], "target": e[1], "type": "modification", "tag": e[2]}

    edges = []
    for e in iso_edges:
       edges.append( (e[0], e[1], {"type": "isotope", "tag": e[2]}) )
    for e in adduct_edges:
        edges.append( (e[0], e[1], {"type": "modification", "tag": e[2]}) )

    G = nx.Graph()
    G.add_edges_from( edges )
    subnetworks = [G.subgraph(c).copy() for c in nx.connected_components(G)]

    return subnetworks, peak_dict, edge_dict


def realign_isotopes(sorted_mz_peak_ids, isotope_search_patterns, mz_tolerance=0.01):
    '''To snap isotopic branch. Assume lowest m/z as M0, and re-align other features against M0.
        Because edges in g can be relationship between any pairs. 
        Re-alignment will get them consistent on grid.
        No redundant features are allowed here, whihc are handled in khipu.clean().

    Parameters
    ----------
    sorted_mz_peak_ids : [(mz, peak_id), ...]; must be unique per m/z. khipu.khipu.clean() takes care of that.
    isotope_search_patterns : [ (1.003355, '13C/12C', (0, 0.8)), (3.010065, '13C/12C*3', (0, 0.8)),..]

    Returns
    -------
    A dictionary of {'M0': F1, '13C/12C*2': F11, ...}
    '''
    M0 = sorted_mz_peak_ids[0]
    _d = {'M0': M0[1]}
    for p in sorted_mz_peak_ids[1:]:
        match = get_isotope_pattern_name(p[0] - M0[0], isotope_search_patterns, mz_tolerance)
        if match == 'Unknown':
            _d["? " + p[1]] = p[1]
            print(p)
        else:
            _d[match] = p[1]

    return _d

def realign_isotopes_reverse(sorted_mz_peak_ids, isotope_search_patterns, mz_tolerance=0.01):
    '''To snap isotopic branch. Assume lowest m/z as M0, and re-align other features against M0.
        Because edges in g can be relationship between any pairs. 
        Re-alignment will get them consistent on grid.
        No redundant features are allowed here, whihc are handled in khipu.clean().

    Parameters
    ----------
    sorted_mz_peak_ids : [(mz, peak_id), ...]; unique per m/z not required, different from realign_isotopes.
    isotope_search_patterns : [ (1.003355, '13C/12C', (0, 0.8)), (3.010065, '13C/12C*3', (0, 0.8)),..]

    Returns
    -------
    A dictionary of {F0: 'M0', F1, '13C/12C*2', ...}
    '''
    M0 = sorted_mz_peak_ids[0]
    _d = {'M0': M0[1]}
    for p in sorted_mz_peak_ids[1:]:
        _d[p[1]] = get_isotope_pattern_name(p[0] - M0[0], isotope_search_patterns, mz_tolerance)
    return _d

def get_isotope_pattern_name(mz, isotope_search_patterns, mz_tolerance=0.01):
    '''Get the isotope with closest m/z match. 
    If error > mz_tolerance, return Unknown, 
    which can happen if the isotope_search_patterns does not cover all possible labled atoms.
    The mz_tolerance needs not to be too precise, as input value was from isotopic_edges.
    Used by realign_isotopes.
    Returns
    -------
    A name in isotope_search_patterns or 'Unknown'.
    '''
    L = sorted([(abs(mz-x[0]), x[1]) for x in isotope_search_patterns])
    if L[0][0] > mz_tolerance:
        print("Warning no match in isotope_pattern: ", mz)
        return 'Unknown'
    else:
        return L[0][1]


def add_data_to_tag(trees, len_limit=20):
    '''
    Append relation note in data to tree.tag, for more informtive showing
    '''
    for tree in trees:
        for node in tree.nodes:
            N = tree.get_node(node)
            N.tag += ' ' + str(N.data)[:len_limit]

    return trees


def export_json_trees(trees, outfile="export_annoTree.tsv"):
    L = [tree.to_json() for tree in trees]
    with open(outfile, 'w', encoding='utf-8') as f:
        json.dump(L, f, ensure_ascii=False, indent=2)

def export_tsv_trees(trees, outfile="export_annoTree.tsv"):
    s = 'Feature_ID\tFeature_tag\troot\troot_tag\trelation\n'
    for tree in trees:
        for node in tree.nodes:
            ND, ROOT = tree.get_node(node), tree.get_node(tree.root)
            s += '\t'.join([ND.identifier, ND.tag, ROOT.identifier, ROOT.tag, ND.data or '']) + '\n'
    with open(outfile, 'w') as O:
        O.write(s)

def is_datatag_in_tree(tree, datatag="13C/12C*6"):
    _in_ = False
    for node in tree.nodes:
        if datatag in tree.get_node(node).tag or datatag in tree.get_node(node).data:
            _in_ = True
            
    return _in_

def find_trees_by_datatag(trees, datatag="13C/12C*6"):
    found = []
    for tree in trees:
        if is_datatag_in_tree(tree, datatag):
            found.append(tree)
            
    return found

def find_trees_by_datatag_list(trees, datatag_list=["13C/12C", 
                                "13C/12C*2", "13C/12C*3", "13C/12C*4", "13C/12C*5", "13C/12C*6",]):
    '''
    Return a list of [tree roots] corresponding to the datatag_list.
    Note 13C/12C is not limited to *1 but all inclusive.
    '''
    return [find_trees_by_datatag(trees, datatag) for datatag in datatag_list]


def __mass_pair_mapping(list1, list2):
    '''Placeholder. Find best matched value_index in list2 for each value in list1.
    Parameters
    ----------
    list1, list2 : each a list of m/z values
    Returns
    -------
    list1_matched : list of index position of best match in list2 for each value in list1
    Notes
    -----
    ref: asari.mass_functions.complete_mass_paired_mapping
    Given the small size of matrices in khipu, np.argmin() is a better solution.
    '''
    pass
