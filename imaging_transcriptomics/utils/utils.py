import os
import warnings
from pathlib import Path

# Third party imports
import numpy as np
with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from netneurotools import freesurfer, stats


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
    Please note that the returned array will have only cortical values (if data corresponding
    to the index of subcortical regions these are automatically discarded by the null permutation
    library). Additionally, we perform permutation on the left hemisphere only (given that the
    right hemisphere of the AHBA is sparsly sampled).

    INPUTS:
    - scan_to_permute: array with the scan values to permute (e.g.: array with average values in some ROIs).
    - n_permutations: number of permutations to do on the scan.
    """
    # Annotation file for the Desikan-Killiany atlas in fs5
    annot_lh = Path(__file__).resolve().parent.parent / \
        "data/fsa5_lh_aparc.annot"
    annot_rh = Path(__file__).resolve().parent.parent / \
        "data/fsa5_rh_aparc.annot"
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
    - user_var: variance explained by the components, selected by the user (deafult= 60%)
    """
    min_explained_variance = user_var
    cumulative_explained_var = np.cumsum(explained_model_variance)
    dim = 1
    while cumulative_explained_var[dim-1] < min_explained_variance:
        dim += 1
    return dim


def pls_alignments(weights, scores):
    """Align PLS components in the desired direction for interpretability.
    """
    weights = -1 * weights
    scores = -1 * scores
    return weights, scores
