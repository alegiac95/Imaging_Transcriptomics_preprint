
==============
How to install
==============

Here is the guide on how to install the *Imaging-Transcriptomics* tool so you can run it from your command line. 
We reccomend using :code:`anaconda` to manage your installation of the :code:`imaging-transcriptomics` package.

To run the script you will need to install two additional packages not listed in :code:`pip`, which are:

* `pyls <https://github.com/rmarkello/pyls/>`_, to perform partial least squares regression
* `spatialnulls <https://markello-spatialnulls.netlify.app/index.html>`_, to perform the permutatiuon of cortical brain areas maintaining the autocorrelation of the brain.

Both this packages are available from Github.

Setting up your environment
---------------------------

A quick way to have the script setup and running is to create a :code:`conda` environment, install the previously mentioned packages and then install the script.

.. warning:: To set up the environment you will need to have `anaconda <https://docs.anaconda.com/anaconda/install/index.html>`_ and `Git <https://git-scm.com/downloads>`_ installed.


An example of how to create the environment is:

.. code-block:: shell

    $ conda create --name imaging_transcriptomics python=3.7 pip
    $ conda activate imaging_transcriptomics
    $ git clone https://github.com/netneurolab/markello_spatialnulls
    $ cd markello_spatialnulls/
    $ conda env update -f environment.yml
    $ cd ..
    $ git clone https://github.com/rmarkello/pyls.git
    $ cd pyls
    $ python setup.py install

The previous list of commands will create a conda environment (named *imaging_transcriptomics*) with the packages needed to run the :code:`imaging_transcriptomics` package.

Installation
------------
To install the package once the previous step is completed run from your termianl the command:

.. code-block:: shell

    $ pip install imaging_transcriptomics

Now you're all set to run your imaging transcriptomics analyses!
Check out our :ref:`getting started <Gettingstarted>` guide and our :ref:`usage <Usage>` to learn how to use the script or the rest of the documentation to learn a bit :ref:`more about the methods <imgtrans>`.


