.. _setup:

=====
Setup
=====

Here we show how to install GPSEA and to prepare your Python environment
for genotype-phenotype association analysis.


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
and install all packages into the active environment.

A specific release (e.g. ``v0.7.0``) can be installed with::

  python3 -m pip install gpsea=0.7.0


Latest release
==============

It is also possible to instal the *latest* release with the latest features and bugfixes.
On top of having Python available, installation of the latest release needs Git to be present as well::

  git clone https://github.com/monarch-initiative/gpsea.git
  
  cd gpsea

  git checkout develop

  python3 -m pip install .


We clone the source code from GPSEA's GitHub repository into a local folder,
then enter the folder, and switch to the `develop` branch that contains the latest features.
Once on `develop`, we use `pip` to use the current folder (`.`) as a package,
and install it into the current Python (virtual) environment.



*****
Tests
*****

Running tests can be done as an optional step after the installation, to ensure GPSEA works correctly.
However, several extra libraries are required to do so, hence we must do one more install::

  python3 -m pip install .[test]

This time, `pip` installs the package with the `test` extras enabled.
Then, running the tests should be as simple as::

  pytest

The Pytest runner will find and run the unit and integration tests, and report the results to the command line.
Only the tests that do *not* require internet access are run.
To run the "online" tests, we add ``--runonline`` option to the command line invocation::

  pytest --runonline

That's all about testing!
