# Virtual Histology

The `virtual_histology` script allows to perform imaging transcriptomics of a neuroimaging file. 

### Usage
The current version can be ran after having checked that you have all necessary dependencies (_see below_, we recommend using `conda` and setup an environment as described in [spatialnulls](https://markello-spatialnulls.netlify.app/setting_up.html)), and downloading the script `virt_hist.py` by calling 

`pyhton virt_hist.py -i /path_to_scan.nii[.gz] -n 1`

or 

`pyhton virt_hist.py -i /path_to_scan.nii[.gz] -v 60`

where `-n` allows to control the number of components for the pls regression while `-v` indicates the percentage of variance to quantify.
The input parameters for the script are:
- -i (--input): path to a neuroimaging file in NIfTI format (.nii,.nii.gz) to analyse. _Please note_ that it is assumed that this file has the same dimensions as the DK atlas used from the script (_1 mm isotropic, matrix size 182x218x182_). If this isn't true the script will raise an error. To resolve this simply reslice your imaging scan to the appropriate dimension.
- -n (--ncomp): number of components to use in the pls regression. _Please note_ that in pls regression the first component is not necessarily the one explaining most variance.
- -v (--variance): minimum explained variance by the pls components. _Please note_ that the amount of variance will be used to calculate the number of pls components (`integer number`) that explains _at least_ the desired variance.
- -o (--output): path where to save the results. if none is given, the results are saved in the same path as the original scan.


__NOTE:__ this initial release (v0.0) is the one tested for development and publication and is _as coded_, while we're currently working on a first release (v1.0) which will have a refarctored/cleaner code and a complete documentation. 
Additionally we are working on an easy installation via `pypi` to allow direct calling from the terminal and in scripts (e.g., bash scripts).

### Requirements
python >= 3.7

[pyls](https://pyls.readthedocs.io/en/latest/)

[spatialnulls](https://markello-spatialnulls.netlify.app/setting_up.html)
