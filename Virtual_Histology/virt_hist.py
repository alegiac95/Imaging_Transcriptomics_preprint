import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    import argparse
    import pandas as pd
    import numpy as np
    import nibabel as nib
    import matplotlib.pyplot as plt
    from pathlib import Path
    from pyls import pls_regression
    from scipy.stats import zscore, norm
    from netneurotools import freesurfer, stats
    from matplotlib import style
    from statsmodels.stats.multitest import multipletests
    from fpdf import FPDF
    from datetime import datetime
    from tqdm import tqdm

# Exceptions defined for the Virtual Histology project.
# The exceptions include most of the possible errors to be dispalyed by the input controls.


class InvalidFormatError(Exception):
    """Exception raised when the format of one of the input files in not correct.

    Attributes:
        errorFile -- file that is not in the correct format.
        message -- optional user overridden error message to display.
    """

    def __init__(self, errorFile, message="The provided file has an invalid format. Please use files in the .nii, .nii.gz format."):
        self.errorFile = errorFile
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.message} The error was caused by the file {self.errorFile}."


class InvalidSizeError(Exception):
    """Exception raised when the size of the images is not correct.

    Attributes:
        errorFile -- file with the wrong size
        size -- size of the input image
        message -- optional user defined error message
    """

    def __init__(self, errorFile, size, message="The provided file has a wrong size."):
        self.errorFile = errorFile
        self.message = message
        self.size = size
        super().__init__(self.message)

    def __str__(self):
        return f"{self.message} The file {self.errorFile} has size: {self.size}"

# CHECKS TO RUN


def check_input_file_exists(file_path):
    """Check if the file given as input to the script exists and if it is in the correct format.
    The possible input formats are NiFTI (.nii) and NiFTI compressed (.nii.gz).
    If either of the two conditions is not met one of 3 possible errors is raised:
    - FileNotFoundError
    - IsADirectoryError
    - InvalidFormat
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError  # Python built-in error
    if file_path.is_dir():
        raise IsADirectoryError  # Implemented in the Errors.py file
    if str().join(file_path.suffixes) not in [".nii.gz", ".nii"]:
        # Implemented in the Errors.py file
        raise InvalidFormatError(file_path)
    return


def check_correct_size(file_path, atlas_path):
    """Check if the dimension of the file is the same as the atlas.
    The used atlas is the Desikan-Killany atlas which is provided in the data directory.
    """
    file_path = Path(file_path)
    atlas_path = Path(atlas_path)
    # Load the scan and atlas using nibabel
    scan = nib.load(file_path)
    atlas = nib.load(atlas_path)
    # Check the sizes of matrixes are the same
    if not atlas.get_fdata().size == scan.get_fdata().size:
        raise InvalidSizeError(file_path, scan.get_fdata().size)
    return

# INPUT PARSING


def get_args():
    """Parse the inputs to run the virtual histology script.
    """
    parser = argparse.ArgumentParser(description="Perform virtual histology analysis on a neuroimaging scan. ",
                                     epilog="Check your results in the specified folder or in the file path of the input scan, if you have not specified an output path.")

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


# PROCESSING FUNCTIONS
def get_average_from_atlas(scan_data, atlas_data, n_regions=41):
    """Extract the average value from the scan in all ROIs from the atlas.
    It is to be noted that only the left hemisphere is used and that the atlas used is the Desikan-Killany
    atlas. This choice hard codes the number of the ROIs that are extraced to 42.

    INPUT:
    - scan_data: array with the scan voxel data (if using nibabel it can be obtained with scan.get_fdata())
    - atlas_data: array with the atlas voxel data (if using nibabel it can be obtained with atlas.get_fdata()).

    Note that the atlas is NOT a binary mask but it is a label file with each value representing a parcellation of the brain,
    either left or right side of the brain.
    """
    data = []
    for i in range(1, n_regions+1):
        data.append(np.mean(scan_data[np.where(atlas_data == i)]))
    data = np.array(data)
    return data


def permute_scan(scan_to_permute, n_permutations=1000):
    """Use null models to permute the scan and get new permuted parcellations.
    The permutations are done with the spatialnulls library available on Github.

    INPUTS:
    - scan_to_permute: array with the scan values to permute (e.g.: array with average values in some ROIs).
    - n_permutations: number of permutations to do on the scan.
    """
    # Annotation file for the Desikan-Killiany atlas in fs5
    annot_lh = Path(__file__).resolve().parent / "data/fsa5_lh_aparc.annot"
    annot_rh = Path(__file__).resolve().parent / "data/fsa5_rh_aparc.annot"
    # Get the parcel centroids of the Desikan-Killiany atlas
    parcel_centroids, parcel_hemi = freesurfer.find_parcel_centroids(
        lhannot=annot_lh, rhannot=annot_rh, version='fsaverage5', surf='sphere', method="surface")
    # Mask the results to have only the left hemisphere
    left_hemi_mask = parcel_hemi == 0
    coords, hemi = parcel_centroids[left_hemi_mask], parcel_hemi[left_hemi_mask]
    # Get the spin samples
    spins = stats.gen_spinsamples(
        coords, hemi, n_rotate=n_permutations, method='vasa', seed=1234)
    y_permute = scan_to_permute[spins]
    return np.array(y_permute)


def get_optimal_components(explained_model_variance, user_var=0.6):
    """Function to return the optimal number of components to explain at least 60% of the explained variance.

    INPUT:
    - explained_model_variance: np.array with the explained variance from a model for each component.
    """
    min_explained_variance = user_var
    cumulative_explained_var = np.cumsum(explained_model_variance)
    dim = 1
    while cumulative_explained_var[dim-1] < min_explained_variance:
        dim += 1
    return dim


def bootstrap_pls(X, Y, Y_permuted, dim, n_permutations=1000):
    """Perform bootstrap on the pls partial least squares regression (PLS).

    INPUTS:
    - X: X data for colleation.
    - Y: original Y data for correlation.
    - Y_permuted: permuted data to use for the correlation.
    - dim: optimal number of components for the PLS regression previously calculated.
    """
    # Initialize the vectors with results
    R_bootstrap = np.zeros(dim)
    P_bootstrap = np.zeros(dim)

    for components in range(1, dim + 1):
        pls_results = pls_regression(
            X, Y, n_components=components, n_perm=0, n_boot=0)
        exp_var = pls_results.get("varexp")
        temp = 100 * exp_var.cumsum(axis=0)
        R_squared = temp[components-1]
        # Perform PLS regression of permuted data
        Rsq = np.zeros(n_permutations)
        for i in tqdm(range(1000), desc=f"Bootstrapping on PLS component {components}", unit="iteration"):
            y_data = Y_permuted[:, i].reshape(41, 1)
            _result = pls_regression(
                X, y_data, n_components=components, n_perm=0, n_boot=0)
            _exp_var = 100 * np.cumsum(_result.get("varexp"))
            Rsq[i] = _exp_var[components-1]
        R_bootstrap[components-1] = R_squared
        P_bootstrap[components-1] = float(len(Rsq[np.nonzero(
            Rsq >= R_squared)])) / n_permutations

    return R_bootstrap, P_bootstrap


def pls_alignments(weights, scores):
    """Align PLS components in the desired direction for interpretability.
    """
    weights = -1 * weights
    scores = -1 * scores
    return weights, scores


# PLOTTING AND REPORTING


class PDF(FPDF):
    def header(self):
        self.rect(10.0, 10.0, 190.0, 280.0)
        self.line(10.0, 50.0, 200.0, 50.0)
        self.set_font("Helvetica", "B", 14)
        self.cell(w=0, h=15.0, align="C",
                  txt="Virtual Histology Analysis Report", ln=True)

    def analysis_info(self, filename, date, filepath):
        self.set_font("Courier", "", 10)
        self.cell(w=100.0, h=8.0, align="L", txt=f"  Scan Name: {filename}")
        self.cell(w=100.0, h=8.0, align="L", txt=f"  Date: {date}", ln=True)
        self.cell(w=100.0, h=10.0, align="L", txt=f"  File Path: {filepath}")

    def pls_regression(self, variance, dim, path_plots):
        self.ln(20)
        self.set_font("Helvetica", "BU", 12)
        self.cell(w=0, h=10.0, align="L", txt="-PLS Regression")
        self.ln(10)
        self.image(Path(path_plots) / "variance.png", w=120.0)
        self.image(Path(path_plots) / "perc_variance.png", w=120.0)


def plot_variance(varexp, dim, sv_path):
    """Generate the plots to use in the pdf.
    The plots are sved as a png to eventually use in presentations and publications.

    INPUTS:
    - varexp: array of variance explained.
    - dim: optimal dimension selected to explain the variance we want.
    - sv_path: path where to save the plots.
    """
    # Plot variance explained --> styling of the plot needed
    save_path = Path(sv_path)
    perc_explained_var = 100 * np.cumsum(varexp)
    plt.plot(range(1, 16), 100 * np.cumsum(varexp),
             marker="o", color="sandybrown")
    plt.plot(dim, perc_explained_var[dim-1], 'o', color="red")
    plt.vlines(
        dim, perc_explained_var[0]-10, perc_explained_var[dim-1], colors="lightgrey", linestyles="dashed")
    plt.hlines(perc_explained_var[dim-1], 0, dim,
               colors="lightgrey", linestyles="dashed")
    plt.title("Cumulative variance explained by PLS components")
    plt.ylabel("Total explained variance (%)")
    plt.xlabel("Number of PLS components")
    plt.xlim(0, 15)
    plt.ylim(perc_explained_var[0]-10, 105)
    plt.grid(True)
    plt.savefig(save_path / "perc_variance.png", dpi=1200)
    plt.close()

    plt.bar(range(1, 16), 100 * varexp, color="sandybrown")
    for index, value in enumerate(varexp):
        plt.text(index + 0.5, 100 * value, "{:.1f}".format(100*value))
    plt.xlabel("PLS components")
    plt.ylabel("Variance (%)")
    plt.savefig(save_path / "variance.png", dpi=1200)
    plt.close()

    return


def reporting(scan_path,  varexp, dim, z, p_val, p_val_corr, pls, output_path=None):
    """Create the reporting for virtual histology.
    The repoting genereated is included in a folder named vh_{filename} and includes:
    - pdf with pls regression components and variance explained by the components.
    - individual csv files with z scores, ROI values, gene ids ans p values for each


    INPUTS:
    - scan_path: path of the scan used for the analysis. This is needed to generate the report folder name.
    - output_path: path where the results are saved. If no path is provided the pathof the scan is used.
    - **kwargs: dictionaries of the values to save to file.
    """

    print("Creating report...")
    if output_path == None:
        output_path = Path(scan_path).parent
    else:
        output_path = Path(output_path)
    scan_path = Path(scan_path)

    scan_name = scan_path.name.split(".")[0]
    save_folder_name = f"vh_{scan_name}"
    matches = list(output_path.glob(f"{save_folder_name}*"))
    if matches:
        nr = len(matches)
        save_folder_name = f"{save_folder_name}_{nr}"
    out_directory = output_path / save_folder_name
    out_directory.mkdir()

    # Create the plots
    plot_variance(varexp, dim, out_directory)

    report_path = out_directory / "report.pdf"
    date = datetime.now().strftime("%d-%m-%Y")
    report = PDF(orientation="P", unit="mm", format="A4")
    report.add_page()
    report.analysis_info(filename=scan_path.name,
                         date=date, filepath=scan_path)
    report.pls_regression(0, 0, path_plots=out_directory)
    report.output(report_path, 'F')

    for i in z.keys():
        data = np.vstack((pls[i].reshape(1, 15633), z[i].reshape(
            1, 15633), p_val[i].reshape(1, 15633), p_val_corr[i].reshape(1, 15633))).T
        data = pd.DataFrame(data, columns=["Gene ID", "Z", "p", "p corrected"])
        data.to_csv(f"{out_directory}/PLS{i}.csv", index=False)


def print_table(var, p):
    print("+-----------+----------------+-------+")
    print("| Component | Cumulative var | p val |")
    print("|-----------|----------------|-------|")
    for i in range(var.size):
        print("|     {}     |      {:.3f}    | {} |".format(
            i+1, var[i], p[i]))
    print("+-----------+----------------+-------+")
    print("")

# MAIN VIRTUAL HISTOLOGY SCRIPT


def virtual_histology(data_path, n_comp=None, var=None, out_directory=None):
    style.use("seaborn-colorblind")
    # LOAD ATLAS FILE
    atlas_path = Path(__file__).parent.absolute() / \
        "data/atlas-desikankilliany_new_1mm_MNI152.nii.gz"
    atlas_data = nib.load(atlas_path).get_fdata()
    # LOAD IMAGE SCAN --> assumes dimensions: 182x218x182 (same as the atlas, 1mm MNI152)
    data_path = Path(data_path)
    img_data = nib.load(data_path).get_fdata()

    # Get average ROI from data
    avg_roi = get_average_from_atlas(img_data, atlas_data).reshape((41, 1))
    pd.DataFrame(avg_roi).to_csv('~/Desktop/ROI_data.csv')
    my_data_y = zscore(avg_roi, axis=0, ddof=1)
    pd.DataFrame(my_data_y).to_csv('~/Desktop/z_data.csv')

    # LOAD EXPRESSION AND GENE DATA
    expression_data_file = Path(__file__).parent / \
        "data/expression_separate_normalization.csv"
    expression_data = pd.read_csv(expression_data_file, sep=',')
    my_data_x = expression_data.iloc[0:41, 2:].to_numpy()
    my_data_x = zscore(my_data_x, ddof=1)
    genes_data_file = Path(__file__).parent / \
        "data/genes_label.txt"
    gene_data = pd.read_fwf(genes_data_file, header=None).to_numpy()

    results = pls_regression(my_data_x, my_data_y,
                             n_components=15, n_perm=0, n_boot=0)
    varexp = results.get("varexp")

    if n_comp == None:
        dim = get_optimal_components(varexp, user_var=float(var)/100)
    else:
        dim = n_comp
    # DATA PERMUTATION AND ROTATION
    my_data_cortical = my_data_y[0:34, 0].reshape(34, 1)
    my_data_subcortical = my_data_y[34:, 0].reshape(7, 1)
    print("Creating scan permutations...")
    cortical_permuted = permute_scan(my_data_cortical).reshape(34, 1000)
    # subcortical_permuted --> random shuffling of the values in the subcortical regions
    subcortical_permuted = np.array(
        [np.random.permutation(my_data_subcortical) for _ in range(1000)]).reshape(7, 1000)
    # merge the two permutated regions
    my_permuted_y = np.zeros((my_data_y.shape[0], 1000))
    my_permuted_y[0:34, :] = cortical_permuted
    my_permuted_y[34:, :] = subcortical_permuted
    # Bootstrap to get p_val, looping on all dimensions from 1 to dim
    R_boot, p_boot = bootstrap_pls(my_data_x, my_data_y, my_permuted_y, dim)

    print_table(R_boot, p_boot)

    # Bootstrap the genes list
    gene_index = np.array(list(range(1, 15633+1)))
    results = pls_regression(my_data_x, my_data_y,
                             n_components=dim, n_boot=0, n_perm=0)
    R1 = np.corrcoef(np.hstack((results.get("x_scores"), avg_roi)),
                     rowvar=False)[0, 1:]
    weights = results.get("x_weights")
    scores = results.get("x_scores")
    for i in range(R1.size):
        if R1[i] < 0:
            weights[:, i], scores[:, i] = pls_alignments(
                weights[:, i], scores[:, i])
    pd.DataFrame(weights).to_csv('~/Desktop/weights.csv')
    pd.DataFrame(scores).to_csv('~/Desktop/scores.csv')
    # Store results in descending order
    # create dictionaries as we don't know a priori how many components the algorithm or the user will select
    pls_w = {}
    x = {}
    pls_id = {}
    gene_id = {}
    pls_wz = {}
    _pls_weights_boot = {}

    for i in range(1, dim+1, 1):
        pls_w.update({i: np.sort(weights[:, i-1], kind='mergesort')[::-1]})
        x.update({i: np.argsort(weights[:, i-1], kind='mergesort')[::-1]})
        pls_id.update({i: gene_data[x[i]]})
        gene_id.update({i: gene_index[x[i]]})
        pls_wz.update({i: zscore(pls_w[i], axis=0, ddof=1)})
        _pls_weights_boot.update(
            {i: np.zeros((weights.shape[0], weights.shape[1], 1000))})
    pd.DataFrame().from_dict(pls_w).to_csv('~/Desktop/pls_weights.csv')

    # x shape : 41x15633
    # x_resampled = np.array([np.random.permutation(my_data_x) for _ in tqdm(range(1000), desc="Permuting gene list")])
    # x_res shape : 1000x41x15633
    # BOOTSTRAPPING ON THE GENE LIST
    for i in tqdm(range(1000), desc="Bootstrapping gene list"):
        myresample = np.random.choice(41, size=41)
        x_data = my_data_x[myresample, :]
        y_data = my_data_y[myresample, :]  # my_permuted_y[:, i].reshape(41, 1)
        _results = pls_regression(
            x_data, y_data, n_components=dim, n_perm=0, n_boot=0)
        _weights = _results.get("x_weights")
        for comp in range(1, dim+1):
            _temp = _weights[:, comp-1]
            _new_w = _temp[x[comp]]
            _corr = np.corrcoef(
                np.hstack((pls_w[comp].reshape((15633, 1)), _new_w.reshape((15633, 1)))), rowvar=False)[0, 1:]
            if _corr < 0:
                _new_w = -1 * _new_w
            _pls_weights_boot[comp][:, comp-1, i] = _new_w

    _pls_std = {}
    _temp = {}
    _z = {}
    _ind = {}
    _pls = {}
    p_val = {}
    p_val_corrected = {}
    for comp in range(1, dim+1):
        _pls_std.update(
            {comp: _pls_weights_boot[comp][:, comp-1, :].std(ddof=1, axis=1)})
        t = pls_w[comp]/_pls_std[comp]
        _temp.update({comp: t})
        _z.update({comp: np.sort(_temp[comp], kind='mergesort')[::-1]})
        _ind.update({comp: np.argsort(_temp[comp], kind='mergesort')[::-1]})
        _pls.update({comp: pls_id[comp][_ind[comp]]})
        p_data = norm.sf(abs(_z[comp]))
        p_val.update({comp: p_data})
        _, corrected, _, _ = multipletests(
            p_data[::-1].reshape(1, 15633), method="fdr_bh", is_sorted=True)
        p_val_corrected.update({comp: corrected})
    pd.DataFrame().from_dict(_temp).to_csv('~/Desktop/temp.csv')
    pd.DataFrame().from_dict(_pls_std).to_csv('~/Desktop/std.csv')
    # Report generation and saving
    reporting(data_path, varexp, dim, z=_z, p_val=p_val,
              p_val_corr=p_val_corrected, pls=_pls, output_path=out_directory)


def main():
    # execute the code
    inputs = get_args()
    data_path = inputs.input
    check_input_file_exists(data_path)
    check_correct_size(data_path, Path(__file__).parent.absolute() /
                       "data/atlas-desikankilliany_new_1mm_MNI152.nii.gz")
    var = inputs.variance
    ncomp = inputs.ncomp
    out_directory = inputs.out
    virtual_histology(data_path, n_comp=ncomp, var=var,
                      out_directory=out_directory)


if __name__ == '__main__':
    main()
