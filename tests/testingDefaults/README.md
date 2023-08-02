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

# Arachnodactyly HP:0001166
robot extract --input hp.obo --method BOT --term HP:0001166 \
  convert --output arachnodactyly.hp.obo

# Focal clonic seizure HP:0002266
robot extract --input hp.obo --method BOT --term HP:0002266 \
  convert --output fcs.hp.obo

# Perimembranous ventricular septal defect HP:0011682
robot extract --input hp.obo --method BOT --term HP:0011682 \
  convert --output pvsd.hp.obo

# Hepatosplenomegaly HP:0001433
robot extract --input hp.obo --method BOT --term HP:0001433 \
  convert --output hepatosplenomegaly.hp.obo

# Tubularization of Bowman capsule HP:0032648
robot extract --input hp.obo --method BOT --term HP:0032648 \
  convert --output bowman.hp.obo

# Intercostal muscle weakness HP:0004878
robot extract --input hp.obo --method BOT --term HP:0004878 \
  convert --output intercostal.hp.obo

# Enuresis nocturna HP:0010677
robot extract --input hp.obo --method BOT --term HP:0010677 \
  convert --output enuresis.hp.obo

# Spasticity HP:0001257
robot extract --input hp.obo --method BOT --term HP:0001257 \
  convert --output spasticity.hp.obo

# Chronic pancreatitis HP:0006280
robot extract --input hp.obo --method BOT --term HP:0006280 \
  convert --output cp.hp.obo

# We use a kind of a hack to include both 
# Clinical modifier HP:0012823
robot extract --input hp.obo --method BOT --term HP:0012823 \
  convert --output cm.bot.hp.obo
robot extract --input hp.obo --method TOP --term HP:0012823 \
  convert --output cm.top.hp.obo
  
# Frequency HP:0040279
robot extract --input hp.obo --method BOT --term HP:0040279 \
  convert --output freq.bot.hp.obo
robot extract --input hp.obo --method TOP --term HP:0040279 \
  convert --output freq.top.hp.obo
  
# Mode of inheritance HP:0000005
robot extract --input hp.obo --method BOT --term HP:0000005 \
  convert --output moi.bot.hp.obo
robot extract --input hp.obo --method TOP --term HP:0000005 \
  convert --output moi.top.hp.obo

# Past medical history HP:0032443
robot extract --input hp.obo --method BOT --term HP:0032443 \
  convert --output pmh.bot.hp.obo
robot extract --input hp.obo --method TOP --term HP:0032443 \
  convert --output pmh.top.hp.obo
  
# Blood group HP:0032223
robot extract --input hp.obo --method BOT --term HP:0032223 \
  convert --output bg.bot.hp.obo
robot extract --input hp.obo --method TOP --term HP:0032223 \
  convert --output bg.top.hp.obo

# Merge into one file
robot merge --input arachnodactyly.hp.obo \
  --input fcs.hp.obo \
  --input pvsd.hp.obo \
  --input hepatosplenomegaly.hp.obo \
  --input bowman.hp.obo \
  --input intercostal.hp.obo \
  --input enuresis.hp.obo \
  --input spasticity.hp.obo \
  --input cp.hp.obo \
  --input cm.bot.hp.obo \
  --input cm.top.hp.obo \
  --input freq.bot.hp.obo \
  --input freq.top.hp.obo \
  --input moi.bot.hp.obo \
  --input moi.top.hp.obo \
  --input pmh.bot.hp.obo \
  --input pmh.top.hp.obo \
  --input bg.bot.hp.obo \
  --input bg.top.hp.obo \
  --output hp.toy.json

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
