import os
import sys
import numpy as np
import nibabel as nib
from pathlib import Path


def normalize_scan(image_path):
    image_path = Path(image_path)
    image_scan = nib.load(image_path)
    file_name = image_path.stem.strip('.nii')
    base_path = image_path.parent
    new_file_path = base_path / f"normalized_{file_name}_1.nii.gz"
    header = image_scan.header
    data = image_scan.get_fdata()
    print(f"min: {data.min()}, max: {data.max()}")
    affine = image_scan.affine
    int_data = data - data.min()
    normalized_data = int_data / np.percentile(int_data, 99)
    print(
        f"NORMALIZED  min: {normalized_data.min()}, max: {normalized_data.max()}")
    new_nifti = nib.Nifti1Image(normalized_data, affine, header=header)
    print(f"{new_file_path}")
    new_nifti.to_filename(new_file_path)


if __name__ == '__main__':
    file_scan = sys.argv[1]
    normalize_scan(file_scan)
