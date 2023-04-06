# gpc
Genotype Phenotype Correlation 


# Set up
create a virtual envirnment. We will then need to install a specific version of setuptools
```
python3 -m venv venv
source venv/bin/activate
pip install "setuptools<58" --upgrade
pip install varcode
```

Then run this from the shell in which the venv has been activate (pysembl gets installed together with varcode)
pyensembl install --release 107 --species homo_sapiens



# Set up dependencies
```
pip install phenopackets
pip install pandas
pip install hpo-toolkit
pip install scipy
pip install statsmodels
```


## Running notebooks

```
python3 -m venv venv
source venv/bin/activate
pip install .
pip install jupyter
python -m ipykernel install --name "venv" --display-name "venv"
cd notebooks
jupyter-notebook
```