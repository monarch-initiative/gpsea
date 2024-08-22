.. _installation:

============
Installation
============

The document describes how to install `gpsea` package into your environment.

Stable release
**************

Installing `gpsea` is easy - we publish releases on `Python Package Index (PyPi) <https://pypi.org/project/gpsea>`_.

Run the following to install the latest *stable* release::

  python3 -m pip install gpsea


Latest release
**************

The *latest* release can be installed by cloning the GitHub repository::

  git clone https://github.com/monarch-initiative/gpsea.git
  cd gpsea

  # Switch to `develop` branch to access the latest features
  git checkout develop

  python3 -m pip install .

The code above will clone the source code from GitHub repository, switch to the `develop` branch with the latest features,
and install the library into the current Python (virtual) environment.


Run tests
^^^^^^^^^

Running tests can be done as an optional step after installation. However, some additional
libraries are required to run tests, hence we must do one more install, this time with `test` option enabled::

  python3 -m pip install .[test]

Then, running the tests is as simple as::

  pytest

This will run the unit and integration tests that do *not* require internet access. To run the "online" tests,
we add ``--runonline`` option to the command line invocation::

  pytest --runonline

That's all about testing!
