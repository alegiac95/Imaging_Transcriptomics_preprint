# Virtual Histology

The `virtual_histology` script allows to perform imaging transcriptomics of a neuroimaging file. 
The current version can be ran after having checked that you have all necessary dependencies (_see below_, we recommend using `conda` and setup an environment as described in [spatialnulls](https://markello-spatialnulls.netlify.app/setting_up.html)) by calling 

`pyhton virt_hist.py -i /path_to_scan.nii[.gz] -n 1`

or 

`pyhton virt_hist.py -i /path_to_scan.nii[.gz] -v 60`

where `-n` allows to control the number of components for the pls regression while `-v` indicates the percentage of variance to quantify.


__NOTE:__ this initial release (v0.0) is the one tested for development and publication and is _as coded_, while we're currently working on a first release (v1.0) which will have a refarctored/cleaner code and a complete documentation. 
Additionally we are working on an easy installation via `pypi` to allow direct calling from the terminal.

### Requirements
python >= 3.7

[pyls](https://pyls.readthedocs.io/en/latest/)

[spatialnulls](https://markello-spatialnulls.netlify.app/setting_up.html)
