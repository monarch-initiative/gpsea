# README

The folder contains JSON files returned by Variant Validator REST API.

## Test `VVMultiCoordinateService`

### `tx_id_PTPN11.json`

Corresponds to `NM_002834.5` in *PTPN11*:

```shell
curl -X 'GET' \
  'https://rest.variantvalidator.org/VariantValidator/tools/gene2transcripts/NM_002834.5' \
  -H 'Content-type:application/json' | python3 -m json.tool > txid-NM_002834.5-PTPN11.json
```

### `tx_id_HBB`

Corresponds to `NM_000518.4` in *HBB*:

```shell
curl -X 'GET' \
  'https://rest.variantvalidator.org/VariantValidator/tools/gene2transcripts/NM_000518.4' \
  -H 'Content-type:application/json' | python3 -m json.tool > txid-NM_000518.4-HBB.json
```

## Test `GeneCoordinateService` functionality of `VVMultiCoordinateService`

### `gene-HBB.json`

Corresponds to all transcripts of *HBB*:

```shell
curl -X 'GET' \
  'https://rest.variantvalidator.org/VariantValidator/tools/gene2transcripts/HGNC:4827' \
  -H 'Content-type:application/json' | python3 -m json.tool > gene-HBB.json
```


### `gene-PTPN11.json`

Corresponds to all transcripts of *PTPN11*:

```shell
curl -X 'GET' \
  'https://rest.variantvalidator.org/VariantValidator/tools/gene2transcripts/HGNC:9644' \
  -H 'Content-type:application/json' | python3 -m json.tool > gene-PTPN11.json
```
