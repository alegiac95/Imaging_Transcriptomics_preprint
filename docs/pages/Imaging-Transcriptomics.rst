
================================
What is imaging transcriptomics?
================================

Imaging transcriptomics is a methodology that allows to identify the spatial autocorrelation between gene expression patterns and some property of brain structure or function as measured by neuroimaging (e.g., MRI, fMRI, PET).

An overview of the methodology can be seen in the figure below.

.. image:: /_static/imaging_transcriptomics.png
    :alt: alternate-text
    :align: center


In brief, average values of the scan are extracted from 41 brain regions as defined by the Desikan-Killiany (DK) atlas. Regional values are then used to perform partial least squares (PLS) regression with gene expression data from the Allen Human Brain Atlas (AHBA) mapped to the DK atlas, in the left hemisphere only. 

As a result of the PLS regression we obtain the ranked genes list according to the spatial alignment with the neuroimaging marker of interest.

.. seealso:: For a more comprehensive dive into the methodology have a look at our paper: 


Allen Human Brain Atlas
-----------------------
The Allen Human Brain Atlas is a 


Desikan-Killiany Atlas
----------------------
The Desikan-Killiany atlas is a parcellation atlas of the human brain, which includes both cortical and subcortical regions.

Partial least squares
---------------------
 


