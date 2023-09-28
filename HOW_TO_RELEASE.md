# How to release

The document describes how to release `genophenocorr` to *PyPi*.

## Release checklist

- update documentation and notebooks to present the latest features 
- create and checkout a release branch
- remove deprecated methods targeted for removal in this version. The `TODO` markers are labeled using 
  the target version (e.g. `TODO[v0.3.0]`)
- bump versions to a release:
  - `src/genophenocorr/__init__.py`
  - `docs/conf.py`
- ensure the CI passes
- deploy to PyPi (described below)
- merge to `main`
- create a GitHub release from the latest `main` commit, including a new tag 
- merge `main` to `develop`
- bump versions to a `dev` version to begin next development iteration

## Deploy to PyPi

### Virtual environment for deployment

As an optional one-time step, consider creating a dedicated virtual environment with the packages required 
for testing, building, and deployment:

```shell
# Create and activate the virtual environment
python3 -m venv build
source build/bin/activate

# Install the build packages
python3 -m pip install build twine
```

### Setup PyPi credentials

As another one-time action, you must create a profile on PyPi and generate an access token. 
The token is used to upload the packages.

First, create an account (e.g. associated with your GitHub account), configure 2FA, store recovery codes, etc.
Then, generate a token to upload the packages. You can generate a token per project or for all projects. 
Store the token into `$HOME/.pypirc` file with `-rw-------` permissions. The file should look like:

```
[pypi]
  username = __token__
  password = <YOUR-TOKEN-HERE>
```

Now we're ready to publish packages!

### Deploy
Run the following to deploy `genophenocorr` to PyPi:  

```bash
# Ensure you're on the release branch
cd genophenocorr

# Build the package
python3 -m build

# Deploy
python3 -m twine upload dist/*

# Clear the built and deployed files
rm -rf build dist
```

The commands will build source distribution and a wheel, and deploy the source distribution and wheel to PyPi.
