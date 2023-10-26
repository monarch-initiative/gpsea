# Variant Validator responses

The folder contains JSON files with example responses from Variant Validator REST API.

## Test `VVHgvsVariantCoordinateFinder`

### `hgvs_MC4R_missense.json`

Corresponds to variant `NM_005912.3:c.253A>G` in *MC4R*:

```shell
curl -X 'GET' \
  'https://rest.variantvalidator.org/VariantValidator/variantvalidator/hg38/NM_005912.3:c.253A>G/mane_select' \
  -H 'Content-type:application/json' | python3 -m json.tool > hgvs_MC4R_missense.json
```

### `hgvs_SURF2_del.json`

Corresponds to variant `NM_017503.5:c.334_336del`:

```shell
curl -X 'GET' \
  'https://rest.variantvalidator.org/VariantValidator/variantvalidator/hg38/NM_017503.5:c.334_336del/mane_select' \
  -H 'Content-type:application/json' | python3 -m json.tool > hgvs_SURF2_del.json
```

### `hgvs_SURF2_dup.json`

Corresponds to variant `NM_017503.5:c.334_336dup`:

```shell
curl -X 'GET' \
  'https://rest.variantvalidator.org/VariantValidator/variantvalidator/hg38/NM_017503.5:c.334_336dup/mane_select' \
  -H 'Content-type:application/json' | python3 -m json.tool > hgvs_SURF2_dup.json
```

### `hgvs_MC4R_dup.json`

Corresponds to variant `NM_005912.3:c.1_3dup` in *MC4R*:

```shell
curl -X 'GET' \
  'https://rest.variantvalidator.org/VariantValidator/variantvalidator/hg38/NM_005912.3:c.1_3dup/mane_select' \
  -H 'Content-type:application/json' | python3 -m json.tool > hgvs_MC4R_dup.json
```

### `hgvs_SURF2_ins.json`

Corresponds to variant `NM_017503.5:c.333_334insAGC` in *SURF2*:

```shell
curl -X 'GET' \
  'https://rest.variantvalidator.org/VariantValidator/variantvalidator/hg38/NM_017503.5:c.333_334insAGC/mane_select' \
  -H 'Content-type:application/json' | python3 -m json.tool > hgvs_SURF2_ins.json
```
