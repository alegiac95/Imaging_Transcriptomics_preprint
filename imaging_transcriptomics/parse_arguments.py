import argparse


def get_args():
    """Parse the inputs to run the virtual histology script.
    """
    parser = argparse.ArgumentParser(description="Perform virtual histology analysis on a neuroimaging scan. ",
                                     epilog="Check your results in the specified folder or in the file path of the input scan, if you have not specified an output path.")

    parser.add_argument("-i", "--input",
                        type=str,
                        help="Input file.",
                        required=True)
    parser.add_argument("-o", "--out",
                        type=str,
                        help="Path where to save the output, if not specified the path of the path of the input scan is used.",
                        required=False)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-n", "--ncomp",
                       type=int,
                       choices=range(1, 16, 1),
                       help="Number of PLS components to use.")
    group.add_argument("-v", "--variance",
                       type=float,
                       choices=range(10, 100),
                       help="""Variance exlplained by the components. 
                       The variance input should be between 10 and 100, and the program will select the number of components that explain a variance closest to the desired (with the explained variance used as a minimum). """)
    return parser.parse_args()
