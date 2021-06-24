import os

# Third party imports
import nibabel as nib
from pathlib import Path

# Custom imports
from .errors import InvalidSizeError, InvalidFormatError


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
