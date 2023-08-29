# Toy HPO ontology

Prepare `hpo.toy.json`, a HPO module that consists of selected terms from different organ systems.

The module contains following descendants of *Phenotypic abnormality* and their ancestors:
- Arachnodactyly HP:0001166
- Focal clonic seizure HP:0002266
- Perimembranous ventricular septal defect HP:0011682
- Hepatosplenomegaly HP:0001433
- Tubularization of Bowman capsule HP:0032648
- Intercostal muscle weakness HP:0004878
- Enuresis nocturna HP:0010677
- Spasticity HP:0001257
- Chronic pancreatitis HP:0006280

On top of *Phenotypic abnormality* descendants, the module contains the *Phenotypic abnormality* siblings 
(e.g. *Clinical modifier*, *Frequency*). 

Prepare the toy JSON by running the following [robot](https://robot.obolibrary.org) commands:

```shell
HPO=https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2023-01-27/hp.obo
module load robot/1.8.3

wget $HPO
```

Extract selected phenotypic features:

```shell
PHENOTYPE_FEATURES=(
'HP:0001166' # Arachnodactyly
'HP:0002266' # Focal clonic seizure
'HP:0011682' # Perimembranous ventricular septal defect
'HP:0001433' # Hepatosplenomegaly
'HP:0032648' # Tubularization of Bowman capsule
'HP:0004878' # Intercostal muscle weakness
'HP:0010677' # Enuresis nocturna
'HP:0001257' # Spasticity
'HP:0006280' # Chronic pancreatitis
)

for term in ${PHENOTYPE_FEATURES[@]}; do
  echo "Processing ${term}"
  robot extract --input hp.obo --method BOT --term ${term} \
  convert --output ${term}.out.hp.obo
done
```

Extract other HPO branches:

```shell
OTHER_BRANCHES=(
'HP:0012823' # Clinical modifier 
'HP:0040279' # Frequency 
'HP:0000005' # Mode of inheritance 
'HP:0032443' # Past medical history 
'HP:0032223' # Blood group 
)

for term in ${OTHER_BRANCHES[@]}; do
  echo "Processing ${term}"
  robot extract --input hp.obo --method BOT --term ${term} \
    convert --output ${term}.bot.out.hp.obo
  robot extract --input hp.obo --method TOP --term ${term} \
    convert --output ${term}.top.out.hp.obo
done
```

Merge into one file:

```shell
INPUTS=""
for obofile in $(ls *.out.hp.obo); do
  INPUTS="--input ${obofile} ${INPUTS}"
done
robot merge ${INPUTS} --output hp.toy.json

rm *.obo
```

# Small HPO annotations file

`phenotype.real-shortlist.hpoa` contains 2 diseases from the HPO annotation file.

# Small HPO ontology

The ontology that contains the terms used in `phenotype.real-shortlist.hpoa` annotation file.

```shell
robot extract --input-iri https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2023-04-05/hp-base.owl \
  -T hp.small.term_ids.txt -o hp.small.owl --method BOT --copy-ontology-annotations true
obographs convert -f json hp.small.owl
rm hp.small.owl
```
