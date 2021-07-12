import setuptools

setuptools.setup(
    name="imaging-transcriptomics",
    version="1.0",
    author="Daniel Martins, MD, PhD; Alessio Giacomel",
    author_email=["daniel.martins@kcl.ac.uk", "alessio.giacomel@kcl.ac.uk"],
    description="Perform imaging-transcriptomics analysis on a neuroimaging scan.",
    # TODO: add long description
    classifiers=["Intended Audience :: Healthcare Industry",
                 "Intended Audience :: Science/Research",
                 "Topic :: Scientific/Engineering :: Image Processing",
                 "Topic :: Scientific/Engineering :: Medical Sciences App", "Development Status :: 4 - Beta",
                 "Programming Language :: Python :: 3",
                 "Programming Language :: Python :: 3.6",
                 "License :: OSI Approved :: MIT License"
                 ],
    keywords="image analysis, neuroimaging, imaging transcriptomics, medical image, research, multimodal imaging",
    python_requires=">=3.6",
    entry_points={"console_scripts": ["imaging-transcriptomics = Imaging_Transcriptomics.imaging_transctiptomics.transcriptomics:main"]})
