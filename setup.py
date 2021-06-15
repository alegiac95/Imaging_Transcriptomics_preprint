import setuptools

setuptools.setup(
    name="Virtual Histology",
    version="0.0",
    author=["Daniel Martins, MD, PhD", "Alessio Giacomel"],
    author_email=["daniel.martins@kcl.ac.uk", "alessio.giacomel@kcl.ac.uk"],
    description="Perform Virtual Histology on the input scan.",
    entry_points={"console_scripts": ["virtual_histology = Virtual_Histology.Virtual_Histology:main"]})
