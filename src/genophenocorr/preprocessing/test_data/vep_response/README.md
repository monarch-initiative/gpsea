# VEP responses

The folder contains JSON files with example VEP responses.

## Test `VepFunctionalAnnotator`

### `missense.json`

Corresponds to variant in `../misssense_test.json`.
```shell
curl 'https://rest.ensembl.org/vep/human/region/16:89279135-89279135/C?LoF=1&canonical=1&domains=1&hgvs=1&mutfunc=1&numbers=1&protein=1&refseq=1&mane=1&transcript_version=1&variant_class=1' \
-H 'Content-type:application/json' | python3 -m json.tool > missense.json
```

### `deletion.json`

Corresponds to variant in `../deletion_test.json`.

```shell
curl 'https://rest.ensembl.org/vep/human/region/16:89284129-89284134/C?LoF=1&canonical=1&domains=1&hgvs=1&mutfunc=1&numbers=1&protein=1&refseq=1&mane=1&transcript_version=1&variant_class=1' \
-H 'Content-type:application/json' | python3 -m json.tool > deletion.json
```


## Test `VepHgvsVariantCoordinateFinder`

### `hgvs_missense.json`

Corresponds to variant `NM_005912.3:c.253A>G` in *MC4R*:

```shell
curl 'https://rest.ensembl.org/vep/human/hgvs/NM_005912.3:c.253A>G?refseq=1' \
-H 'Content-type:application/json' | python3 -m json.tool > hgvs_MC4R_missense.json
```

### `hgvs_SURF2_del.json`

Corresponds to variant `NM_017503.5:c.334_336del`:

```shell
curl 'https://rest.ensembl.org/vep/human/hgvs/NM_017503.5:c.334_336del?refseq=1' \
-H 'Content-type:application/json' | python3 -m json.tool > hgvs_SURF2_del.json
```

### `hgvs_SURF2_dup.json`

Corresponds to variant `NM_017503.5:c.334_336dup`:

```shell
curl 'https://rest.ensembl.org/vep/human/hgvs/NM_017503.5:c.334_336dup?refseq=1' \
-H 'Content-type:application/json' | python3 -m json.tool > hgvs_SURF2_dup.json
```

### `hgvs_MC4R_dup.json`

Corresponds to variant `NM_005912.3:c.1_3dup` in *MC4R*:

```shell
curl 'https://rest.ensembl.org/vep/human/hgvs/NM_005912.3:c.1_3dup?refseq=1' \
-H 'Content-type:application/json' | python3 -m json.tool > hgvs_MC4R_dup.json
```

### `hgvs_SURF2_ins.json`

Corresponds to variant `NM_017503.5:c.333_334insAGC` in *SURF2*:

```shell
curl 'https://rest.ensembl.org/vep/human/hgvs/NM_017503.5:c.333_334insAGC?refseq=1' \
-H 'Content-type:application/json' | python3 -m json.tool > hgvs_SURF2_ins.json
```
