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



# Set up phenopackets
```
pip install phenopackets
pip install pandas
```
