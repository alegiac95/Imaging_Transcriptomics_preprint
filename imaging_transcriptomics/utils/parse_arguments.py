import argparse


def get_args():
    """Parse the inputs to run the virtual histology script.
    """
    DESCRIPTION = """Perform virtual histology analysis on a neuroimaging scan. """
    EPILOG = """Check your results in the specified folder or in the file path of the input scan, if you have not specified an output path.
            If you used this software in your research please cite:
            - Imaging transcriptomics: Convergent cellular, transcriptomic, and molecular neuroimaging signatures in the healthy adult human brain
                Daniel Martins, Alessio Giacomel, Steven CR Williams, Federico E Turkheimer, Ottavia Dipasquale, Mattia Veronese, PET templates working group
                bioRxiv 2021.06.18.448872; doi: https://doi.org/10.1101/2021.06.18.448872
            """
    parser = argparse.ArgumentParser(description=DESCRIPTION,
                                     epilog=EPILOG)

    parser.add_argument("-i", "--input",
                        type=str,
                        help="""Input imaging file in NiFTI format (.nii, .nii.gz).
                        The input file is expected to have the same matrix size as the atlas used (182x218x182), 
                        if the input image has different matrix size this can be resliced to match the resolution of the MNI152 template provided with FSL.
                        """,
                        required=True)
    parser.add_argument("-o", "--out",
                        type=str,
                        help="Path where to save the output, if not specified the path of the path of the input scan is used.",
                        required=False)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-n", "--ncomp",
                       type=int,
                       help="Number of PLS components to use. The number of components has to be between 1 and 15.")
    group.add_argument("-v", "--variance",
                       type=int,
                       help="""Variance exlplained by the components.
                       The variance input should be between 10 and 100, and the program will select the number of components that explain a variance closest to the desired (with the explained variance used as a minimum). """)
    return parser.parse_args()