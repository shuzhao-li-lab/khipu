'''
Constructing empirical compounds de novo based on khipu package.
The ion patterns here are constructed differently from formula based calculations,
because formulae are not known here.
Search functions are based on mass2chem.

empCpds (or epds) can be initiated by signatures on isotopic relationships or adduct relationships.
For low-intensity peaks, their isotopic counterparts may not be detectable.
For assigned peaks/features, will calculate selectivity/score later too.
'''

# will add coelution_function
# from mass2chem.search import is_coeluted_by_overlap, is_coeluted_by_distance

from .extended import *

class epdsConstructor:
    '''Wrapper class to organize a list of peaks/features into a list of empirical compounds.

    To-dos:
        add support of user input formats where rtime isn't precise or unavailable.
        add options of coelution_function (see mass2chem.epdsConstructor )
    '''
    def __init__(self, peak_list, mode='pos'):
        '''
        Parameters
        ----------
        peak_list : [{'parent_masstrace_id': 1670, 'mz': 133.09702315984987, 'rtime': 654, 'height': 14388.0, 'id': 555}, ...]
        mz_tolerance_ppm: ppm tolerance in examining m/z patterns.
        '''
        self.mode = mode
        self.peak_list = peak_list

    def peaks_to_epdDict(self, 
                        isotope_search_patterns,
                        adduct_search_patterns,
                        extended_adducts,
                        mz_tolerance_ppm,
                        rt_tolerance=2,
        ):
        '''
        Parameters
        ----------
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
            adduct_search_patterns is dependent on ionization, but the option is left open for other functions.

        mz_tolerance_ppm : ppm tolerance in examining m/z patterns.
        rt_tolerance : tolerance threshold for deviation in retetion time, arbitrary unit depending on input data.
                Default intended as 2 seconds.

        Returns
        -------
        epdDict : A dictionary of empCpds (empirical compounds) indexed by IDs ('interim_id').
                Not including singletons.
        '''
        subnetworks, peak_dict, _ = peaks_to_networks(self.peak_list,
                    isotope_search_patterns,
                    adduct_search_patterns,
                    mz_tolerance_ppm,
                    rt_tolerance
        )
        WV = Weavor(peak_dict, isotope_search_patterns=isotope_search_patterns, 
                    adduct_search_patterns=adduct_search_patterns, 
                    mz_tolerance_ppm=mz_tolerance_ppm, mode=self.mode)
        print("\n")
        print("Initial khipu search grid: ")
        print(WV.mzgrid)
        print("\n")

        khipu_list = graphs_to_khipu_list(
            subnetworks, WV, mz_tolerance_ppm=mz_tolerance_ppm,
            )
        khipu_list, all_assigned_peaks = extend_khipu_list(khipu_list, peak_dict, 
                        extended_adducts, mz_tolerance_ppm=mz_tolerance_ppm,
                        rt_tolerance=rt_tolerance)

        print("\n\n ~~~~~~ Got %d khipus, with %d features ~~~~~~~ \n\n" 
                    %(len(khipu_list), len(all_assigned_peaks)))
        empCpds = export_empCpd_khipu_list(khipu_list)
        epdDict = {}
        for E in empCpds:
            epdDict[E['interim_id']] = E
        return epdDict
