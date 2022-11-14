import networkx as nx

from mass2chem.search import build_centurion_tree, find_all_matches_centurion_indexed_list


def make_peak_tag(peak):
    '''
    peak format: {'id': 'F1', 'mz': 60.0808, 'rtime': 117.7, ...,
    'intensities': [250346.0], 'representative_intensity': 250346.0}
    '''
    return str(round(peak['mz'], 4)) + '@' + str(round(peak['rtime'], 1))

def make_edge_tag(edge):
    '''
    Concatenate str edges by underscore
    '''
    if edge[0] > edge[1]:
        return edge[1] + '_' + edge[0]
    else:
        edge[0] + '_' + edge[1]

def make_peak_dict(peak_list):
    '''
    Same as search.build_peak_id_dict but uses 'id' not 'id_number'.
    '''
    peak_dict = {}
    for p in peak_list:
        peak_dict[p['id']] = p
    return peak_dict

def rt_matched_by_tolerance(P1, P2, rt_tolerance):
    return abs(P2['rtime']-P1['rtime']) < rt_tolerance

def rt_compared_by_values(P1, P2, rt_tolerance=None):
    return P2['rtime'] > P1['rtime']

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
    seed_search_patterns: a list in the format of [(mz difference, notion, (ratio low limit, ratio high limit)), ..]
            This can be obtained through search.isotopic_patterns. The ratios are optional, because 
            1) naturally occuring constrains are based on chemical formula;
            2) rules are different when isotope tracers are introduced to the experiments.
            But it's important to have a comprehensive list here for isotope tracing experiments.
    isotope_rt_tolerance: tolerance threshold for deviation in retetion time, arbitrary unit depending on input data.
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
    seed_search_patterns: a list in the format of [(mz difference, notion, (ratio low limit, ratio high limit)), ..]
            This can be obtained through search.isotopic_patterns. The ratios are optional, because 
            1) naturally occuring constrains are based on chemical formula;
            2) rules are different when isotope tracers are introduced to the experiments.
            But it's important to have a comprehensive list here for isotope tracing experiments.
    isotope_rt_tolerance: tolerance threshold for deviation in retetion time, arbitrary unit depending on input data.
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
    '''

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
    for e in iso_edges + adduct_edges:
        edges.append(e[:2])
    G = nx.Graph()
    G.add_edges_from( edges )
    subnetworks = [G.subgraph(c).copy() for c in nx.connected_components(G)]

    return subnetworks, peak_dict, edge_dict








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




