# ion patterns in high-resolution mass spectrometry
#
# m/z differences corresponding to adducts and isotopes
# Not using charge in tables, because all m/z diff is relative to other charged ions
#
# In-source fragments, i.e. neutral loss, are described in 
# Chi, Y., Mitchell, J.M., Zheng, S., Thapa, M. and Li, S., 2026. Assessing the metabolomics “dark matter” by a detectable khipu model. Metabolomics, 22(3), p.79.


PROTON = 1.00727646677
electron = 0.000549

# These are single-charged values. 
# Double and Tripple charged ions are computed after the initial round, on the fly
# avoid confusing adducts in initial search, e.g. H, H2O

# Long form for isotope labelled data
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

# regular data. Ratios considering measurement errors
isotope_search_patterns_pos = [ (1.003355, '13C/12C', (0, 0.5)),
                            (2.00671, '13C/12C*2', (0, 0.2)),
                            (3.010065, '13C/12C*3', (0, 0.1)),
                            # (3.9948, '44Ca/40Ca', (0, 0.1)), # 2%
                            (1.9970, '37Cl/35Cl', (0.1, 0.8)), # 24.24%
                            ]

isotope_search_patterns_neg = [ (1.003355, '13C/12C', (0, 0.5)),
                            (2.00671, '13C/12C*2', (0, 0.2)),
                            (3.010065, '13C/12C*3', (0, 0.1)),
                            (1.9970, '37Cl/35Cl', (0.1, 0.8)), # 24.24%
                            (1.9958, '32S/34S', (0, 0.1)), # 4%
                            ]

adduct_search_patterns = adduct_search_patterns_pos = [  # initial patterns are relative to M+H+
                            (21.9820, 'Na/H'),
                            (41.026549, 'ACN'),     # Acetonitrile
                            (67.987424, 'NaCOOH'),
                            (37.955882, 'K/H'),
                            (32.026215, 'CH3OH')
                            ]

adduct_search_patterns_neg = [ (35.9767, 'HCl'), 
                            (46.00548, 'HCOOH'),
                            (67.987424, 'NaCOOH'),
                            (82.0030, 'NaCH2COOH'),
                            (21.9820, 'Na/H'), 
                            (41.026549, 'ACN'),
                            (37.955882, 'K/H'),
                            ]

# added ISFs from Chi et al, 2026  Metabolomics, 22(3), p.79
extended_adducts = [
                            (1.0078, 'H'),
                            (17.02655, 'NH3'),
                            (18.0106, 'H2O'),      # easy to confuse with bio reactions
                            (18.033823, 'NH4'),
                            (27.01089904, 'HCN'),
                            (27.99492, 'CO'),
                            (32.026215, 'CH3OH'),
                            (35.9767, 'HCl'),
                            (37.94694, 'Ca/H2'),
                            (43.96389, 'Na2/H2'),
                            (46.00548, 'CO2H2'),
                            (67.987424, 'NaCOOH'),
                            (83.961361, 'KCOOH'),
                            (97.96737927, 'H2SO4'),
                            (97.97689507, 'H3PO4'),
                            # ISFs, i.e., neutral loss
                            (-135.9748, '-2C-4O-2H-2Na'),
                             (-88.0522, '-4C-8H-2O'),
                             (-82.0029, '-2C-H-3F'),
                             (-67.9874, '-C-2O-Na-H'),
                             (-60.021, '2C4H2O'),
                             (-46.0054, '-C-2H-2O'),
                             (-44.0261, '-2C-4H-O'),
                             (-43.9898, '-C-2O'),
                             (-42.0468, '-3C-6H'),
                             (-42.0104, '-2C-2H-O'),
                             (-41.0263, '-2C-3H-N'),
                             (-35.037, '-5H-N-O'),
                             (-30.0468, '-2C-6H'),
                             (-30.0105, '-C-2H-O'),
                             (-28.0312, '-2C-4H'),
                             (-27.9948, '-C-O'),
                             (-26.0156, '-2C-2H'),
                             (-23.9999, '-2C'),
                             (-18.0104, '-2H-O'),
                             (-17.0264, '-N-3H'),
                             (-16.0312, '-C-4H'),
                             (-15.9948, '-O'),
                             (-14.0156, '-C-2H'),
                             (-13.9792, '2H-O'),
                             (-12.0, '-C'),
                             (-2.0155, '-2H')
]

