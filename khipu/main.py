import argparse
from khipu import __version__
from .extended import *


def main():
    '''
    khipu, annotating metabolomics features to empCpds
    '''
    print("\n\n~~~~~~~ Hello from Khipu (%s) ~~~~~~~~~\n" %__version__)

    parser = argparse.ArgumentParser(description='khipu, annotating metabolomics features to empCpds')
    parser.add_argument('-v', '--version', action='version', version=__version__, 
            help='print version and exit')
    parser.add_argument('-m', '--mode', default='pos', 
            help='mode of ionization, pos or neg')
    parser.add_argument('--ppm', default=5, type=int, 
            help='mass precision in ppm (part per million), same as mz_tolerance_ppm')
    parser.add_argument('--rtol', default=2, type=float, 
            help='tolerance of retention time match, arbitrary unit dependent on preprocessing tool')
    parser.add_argument('-i', '--input', 
            help='input file as feature table')
    parser.add_argument('-o', '--output', 
            help='prefix of output files')

    args = parser.parse_args()
    khipu_annotate(args)


#
# -----------------------------------------------------------------------------
#
if __name__ == '__main__':
    
    main()
