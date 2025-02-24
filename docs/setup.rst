.. _setup:

#####
Setup
#####

Here we show how to install GPSEA and prepare your Python environment
for genotype-phenotype association analyses.


.. contents:: Table of Contents
  :depth: 1
  :local:


*************
Install GPSEA
*************

We provide GPSEA in two versions: *stable* with stabilized APIs
and *latest* with bleeding edge features and bugfixes.


Stable release
==============

Most users should install the *stable* version. The installation is easy - we publish the releases
on `Python Package Index (PyPi) <https://pypi.org/project/gpsea>`_.
Therefore, we can use `pip`, the Python package manager, to install latest GPSEA release 
into the active Python (virtual) environment::

  python3 -m pip install gpsea


`pip` will fetch GPSEA from PyPi along with all its dependencies (e.g. SciPy, Matplotlib)
and install the packages into the environment.

A specific release (e.g. ``v0.7.0``) can be installed with::

  python3 -m pip install gpsea=0.7.0


Latest release
==============

It is also possible to instal the *latest* release with the latest features and bugfixes.
On top of having Python available, installation of the latest release needs
`Git <https://git-scm.com/>`_ to be installed as well::

  git clone https://github.com/P2GX/gpsea.git
  
  cd gpsea

  git checkout develop

  python3 -m pip install .


The commands clone the source code from GPSEA's GitHub repository into a local folder,
then we enter the folder and switch to the `develop` branch that contains the latest features.
Once on `develop`, we use ``pip`` to use the current folder (``.``) as a package,
and install it into the current Python (virtual) environment.


*********
Run tests
*********

Running tests can be done as an optional step after the installation, to ensure GPSEA works correctly.
However, several extra libraries are required to do so, hence we must install GPSEA with ``test`` extras::

  python3 -m pip install .[test]

.. note::

  Most users do *not* need to run the tests while installing GPSEA,
  as the tests are run within the contiunuous integration (CI) pipeline.

With the ``test`` extras installed, Then, running the tests is as simple as::

  pytest

The Pytest runner will run >400 unit and integration tests, and report the results to the command line.

Only the tests that do *not* require internet access are run.
To run the "online" tests, we add ``--runonline`` option to the command line invocation::

  pytest --runonline

That's all about testing! Now let's move on to the tutorial,
where we show an example genotype-phenotype association analysis.
