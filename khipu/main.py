'''
Placeholder

'''

import argparse
from khipu import __version__
from .utils import khipu_annotate

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
    parser.add_argument('-i', '--input', 
            help='input directory of mzML files to process, or a single file to analyze')
    parser.add_argument('-o', '--output', 
            help='output directory')
    parser.add_argument('-j', '--project', 
            help='project name')

    args = parser.parse_args()

    khipu_annotate(args)


#
# -----------------------------------------------------------------------------
#
if __name__ == '__main__':
    
    main()