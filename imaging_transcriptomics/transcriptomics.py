# Standard library imports
from pathlib import Path

# Third party imports
import pandas as pd
import numpy as np
import nibabel as nib
from pyls import pls_regression
from scipy.stats import zscore, norm
from matplotlib import style
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

# Custom imports
from . import utils

# PROCESSING FUNCTIONS


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
# MAIN IMAGING TRANSCRIPTOMICS SCRIPT


def imaging_transcriptomics(data_path, n_comp=None, var=None, out_directory=None):
    style.use("seaborn-colorblind")
    # LOAD ATLAS FILE
    atlas_path = Path(__file__).parent.absolute() / \
        "data/atlas-desikankilliany_new_1mm_MNI152.nii.gz"
    atlas_data = nib.load(atlas_path).get_fdata()
    # LOAD IMAGE SCAN --> assumes dimensions: 182x218x182 (same as the atlas, 1mm MNI152)
    data_path = Path(data_path)
    img_data = nib.load(data_path).get_fdata()

    # Get average ROI from data
    avg_roi = utils.get_average_from_atlas(
        img_data, atlas_data).reshape((41, 1))
    my_data_y = zscore(avg_roi, axis=0, ddof=1)

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

    # Select the numebr of components if not specified by the user
    if n_comp == None:
        dim = utils.get_optimal_components(varexp, user_var=float(var)/100)
    else:
        dim = n_comp

    # DATA PERMUTATION AND ROTATION
    my_data_cortical = my_data_y[0:34, 0].reshape(34, 1)
    my_data_subcortical = my_data_y[34:, 0].reshape(7, 1)
    print("Creating scan permutations...")

    # Cortical permuted --> permutation keeping in consideration spatial
    # autocorrelation of the cortical regions
    cortical_permuted = utils.permute_scan(my_data_cortical).reshape(34, 1000)

    # Subcortical_permuted --> random shuffling of the values in the subcortical regions
    subcortical_permuted = np.array(
        [np.random.permutation(my_data_subcortical) for _ in range(1000)]).reshape(7, 1000)

    # Merge permuted cortical regions and shuffled subcortical into the semi-random matrix with the data
    # the matrix is not fully random as there is still a subdivision of the corticla and subcortical data.
    my_permuted_y = np.zeros((my_data_y.shape[0], 1000))
    my_permuted_y[0:34, :] = cortical_permuted
    my_permuted_y[34:, :] = subcortical_permuted

    # Bootstrap to get p_val, looping on all dimensions from 1 to dim
    R_boot, p_boot = bootstrap_pls(my_data_x, my_data_y, my_permuted_y, dim)

    utils.print_table(R_boot, p_boot)

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
            weights[:, i], scores[:, i] = utils.pls_alignments(
                weights[:, i], scores[:, i])

    # Create dictionaries as we don't know a priori how many components the algorithm or the user will select
    # the dictionaries allow to dynamically increase the dimension
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
    # Report generation and saving
    utils.reporting(data_path, varexp, dim, z=_z, p_val=p_val,
                    p_val_corr=p_val_corrected, pls=_pls, output_path=out_directory)


def main():
    # Get user inputs and check them
    inputs = utils.get_args()
    data_path = inputs.input
    utils.check_input_file_exists(data_path)
    utils.check_correct_size(data_path, Path(__file__).parent.absolute() /
                             "data/atlas-desikankilliany_new_1mm_MNI152.nii.gz")
    var = inputs.variance
    ncomp = inputs.ncomp
    out_directory = inputs.out
    # Run the analysis
    imaging_transcriptomics(data_path, n_comp=ncomp, var=var,
                            out_directory=out_directory)


if __name__ == '__main__':
    main()
