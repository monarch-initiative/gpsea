# VEP responses

The folder contains JSON files with example VEP responses.

## `missense.json`

Corresponds to variant in `../misssense_test.json`.
```shell
curl 'https://rest.ensembl.org/vep/human/region/16:89279135-89279135/C?LoF=1&canonical=1&domains=1&hgvs=1&mutfunc=1&numbers=1&protein=1&refseq=1&mane=1&transcript_version=1&variant_class=1' \
-H 'Content-type:application/json' | python3 -m json.tool > missense.json
```

## `deletion.json`

Corresponds to variant in `../deletion_test.json`.

```shell
curl 'https://rest.ensembl.org/vep/human/region/16:89284129-89284134/C?LoF=1&canonical=1&domains=1&hgvs=1&mutfunc=1&numbers=1&protein=1&refseq=1&mane=1&transcript_version=1&variant_class=1' \
-H 'Content-type:application/json' | python3 -m json.tool > deletion.json
```