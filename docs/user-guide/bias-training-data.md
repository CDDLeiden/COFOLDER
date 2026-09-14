# Providing bias training data

COFOLDER's database-backed bias workflow takes searchable source databases, not
query-specific reference CSV files:

```text
--bias_training_data_protein_path ~/.cofolder/data/bias/protein
--bias_training_data_ligand_path  ~/.cofolder/data/bias/ligand
```

These are the defaults, so a normal run can simply use `--assess_bias` after the
databases have been prepared. Database preparation uses network access; assessment
and `--preflight_only` are offline.

## Prepare the standard PDB sources

Install COFOLDER's analysis dependencies and MMseqs2, then run:

```bash
cofolder-tools fetch-bias-training-data \
  --output_root ~/.cofolder/data/bias
```

This downloads the wwPDB Chemical Component Dictionary, asks MMseqs2 to prepare its
`PDB` target database, and snapshots PDB entry/non-polymer metadata from the RCSB
Search and Data APIs. This can require substantial time, network traffic, and disk
space. The command writes complete bundles only after their contents validate:

```text
~/.cofolder/data/bias/
├── protein/
│   ├── manifest.json
│   ├── pdb_release_metadata.csv.gz
│   ├── pdb_sequences.csv.gz
│   └── mmseqs/
│       └── db*
└── ligand/
    ├── manifest.json
    └── pdb_ligands.csv.gz
```

The protein metadata maps PDB identifiers to initial release dates. The checksummed
sequence index is exported from the same MMseqs target and lets plot backfilling
remain offline. The ligand table
contains `pdb_id`, `release_date`, `ligand_id`, and canonical `smiles`. Each manifest
records the source snapshot, schema, checksums, record counts, and tool versions.

Validate an existing download without changing it:

```bash
cofolder-tools fetch-bias-training-data \
  --output_root ~/.cofolder/data/bias \
  --validate_only
```

RCSB entries are retrieved in bounded concurrent batches with retry limits, and each
completed entry is appended to a checkpoint under the output root. An interrupted
preparation therefore retains both the previous complete bundle and resumable fetch
progress. Rerun the command to resume/reuse existing inputs, or pass `--overwrite`
to replace complete bundles.

## Use the sources

The default locations need no path arguments:

```bash
cofolder validate -s system.yaml -o options.yaml --assess_bias
```

Override either database explicitly:

```bash
cofolder validate -s system.yaml -o options.yaml --assess_bias \
  --bias_training_data_protein_path /data/bias/protein \
  --bias_training_data_ligand_path /data/bias/ligand
```

The same arguments work with `bias`, `screen`, and `oracle`. Python recipes expose
the same names:

```python
from cofolder.recipes.validate import Validate

workflow = Validate(
    wrk_dir="runs/example",
    system_path="system.yaml",
    options_path="options.yaml",
    assess_bias=True,
    bias_training_data_protein_path="/data/bias/protein",
    bias_training_data_ligand_path="/data/bias/ligand",
)
workflow.run()
```

Run `--preflight_only` first to validate manifests, selected chains, cutoff, MSA
needs, runner availability, and execution counts without searching the databases or
starting inference.

## Generated per-query references

For each system, COFOLDER searches every distinct selected protein sequence and
ligand. It writes the derived evidence beneath the run rather than modifying the
source bundles:

```text
<wrk_dir>/results/bias_train/
├── reference_manifest.json
├── protein_references.csv
├── ligand_references.csv
├── ligand_references_<CHAIN>.csv
├── bias_training_data.csv
└── reference_landscape_summary.csv
```

`reference_manifest.json` binds the source-database fingerprints to query hashes,
selected chains, thresholds, release policy, and builder schema. Raw searches use a
content-addressed cache at `~/.cofolder/cache/bias_queries` by default. Override it
with `--bias_query_cache_path`; each run manifest records component hit/miss/rebuild
events while every output row retains `query_chain_id`.

Inspect or recoverably invalidate the cache with:

```bash
cofolder-tools bias-query-cache inspect
cofolder-tools bias-query-cache invalidate --kind all
```

Invalidation moves entries to a timestamped sibling directory instead of deleting
an arbitrary path.

One pair of source bundles supports many systems. Distinct protein sequences receive
independent MMseqs searches and distinct ligands receive independent ECFP/Tanimoto
searches. Use `--bias_chains A B` to restrict which query entities are assessed.

## Customize or expand the sources

Prepare bundles outside the default location:

```bash
cofolder-tools fetch-bias-training-data \
  --output_root /data/custom-bias \
  --mmseqs_db_name PDB \
  --bias_training_data_protein_path /data/custom-bias/protein \
  --bias_training_data_ligand_path /data/custom-bias/ligand
```

A custom protein bundle must contain an MMseqs target whose identifiers resolve to
PDB IDs, plus release metadata conforming to the standard bundle. A UniRef-only
database cannot implement PDB release-cutoff semantics unless it also supplies a
documented PDB-ID/release-date mapping. A custom ligand bundle may add occurrences,
but every row must preserve `pdb_id`, `release_date`, `ligand_id`, and `smiles`.

Package an already-built PDB-compatible MMseqs target and its mapping with
`--mmseqs_output_db` and `--protein_release_metadata`. Package a maintained ligand
occurrence table with `--ligand_occurrence_table`. For example:

```bash
cofolder-tools fetch-bias-training-data \
  --output_root /data/custom-bias \
  --mmseqs_output_db /imports/expanded-pdb/db \
  --protein_release_metadata /imports/expanded-pdb/releases.csv.gz \
  --ligand_occurrence_table /imports/pdb-ligands-expanded.csv.gz \
  --overwrite
```

Preparation exports the offline PDB sequence index from the MMseqs target. For an
unusual imported database, supply an existing `pdb_id,sequence` CSV with
`--protein_sequence_table`.

The input tables are schema-checked while the staged bundles are built. The tool
then regenerates manifests and checksums and publishes both bundles atomically;
`--overwrite` never applies to workflow commands. Run `--validate_only` afterward.
Bundle manifests use relative file names, so a whole bundle directory can be moved
or copied without rewriting it.

## Current cutoff policy

The default is:

```text
--bias_release_cutoff 2023-06-01
```

Protein and ligand evidence both use the strict rule
`release_date < 2023-06-01`; an entry released on 2023-06-01 is excluded. Override
the boundary when required:

```bash
cofolder validate -s system.yaml -o options.yaml --assess_bias \
  --bias_release_cutoff 2024-01-01
```

The current protein method searches the MMseqs `PDB` source and records MMseqs
`pident`. The ligand method starts from wwPDB CCD chemistry and PDB occurrence data
from the official [RCSB APIs](https://www.rcsb.org/docs/programmatic-access/web-apis-overview).
The repository does not establish that 2023-06-01 equals every supported backend's
actual training cutoff. Treat it as a configurable reference-overlap boundary, not
as a guaranteed model-training boundary.

Use `--bias_release_cutoff whole` to include every dated entry in the prepared
snapshot. This is a snapshot policy, not a guessed future date and not a statement
about a backend's training cutoff. Exact source maximum dates and snapshot
fingerprints are written to each reference manifest.

## Prepare custom crystal-complex supplements

Copy the packaged examples and edit `custom_bias_complexes.yaml`. Each YAML entry
has a stable `complex_id`, a PDB/mmCIF path, explicit protein chains, and explicit
ligand residue selectors. Ligand SMILES may be supplied directly or resolved from
mmCIF/CCD descriptors; COFOLDER never guesses bonds from coordinates.

```bash
cofolder-tools prepare-bias-custom-complexes custom_bias_complexes.yaml \
  --output_root custom-bias
```

Use the resulting paired bundle with `--custom_bias_reference_path custom-bias`.
Custom structures are supplements and are not filtered by the public release date.
Their manifest fingerprint and `complex_id` pairing provenance are retained.

### Assess PDB references before 2021-09-01 plus custom crystals

The following standalone assessment combines the prepared public PDB snapshot with
the custom crystal bundle:

```bash
cofolder bias \
  --system_path system.yaml \
  --wrk_dir runs/bias-pdb-2021-plus-custom \
  --bias_training_data_protein_path ~/.cofolder/data/bias/protein \
  --bias_training_data_ligand_path ~/.cofolder/data/bias/ligand \
  --custom_bias_reference_path custom-bias \
  --bias_release_cutoff 2021-09-01 \
  --bias_protein_similarity_threshold 0.25 \
  --bias_ligand_similarity_threshold 0.35
```

Public protein and ligand evidence follows the strict rule
`release_date < 2021-09-01`; an entry released on 2021-09-01 is excluded. Custom
complexes supplement that public slice and are not date-filtered. The standalone
`bias` command performs reference-overlap assessment without running a prediction
backend. Add `--preflight_only` to validate the request without searching or writing
outputs.

## Run the dated/whole comparison matrix

`bias_matrix.yaml` runs the four S5.1 comparisons under separate directories:
strictly before 2023-06-01, the whole prepared snapshot, strictly before
2023-06-01 plus the prepared custom supplement, and the custom supplement alone.

```bash
cofolder-tools run-bias-matrix bias_matrix.yaml
```

The installed runner invokes the public `Bias` or `Validate` recipe and writes a
top-level `matrix_manifest.json`; it does not implement a separate scoring path.
The date is an exact comparison cutoff, not a claim about a backend's training date.

When a ligand occurrence belongs to a PDB absent from the selected protein view,
COFOLDER backfills its protein axis from the bundle's offline sequence index using
MMseqs2. Protein and ligand boundaries are independently configurable with
`--bias_protein_similarity_threshold` and
`--bias_ligand_similarity_threshold`, both on the normalized 0–1 plotting scale.
Their defaults are `0.25` (25% MMseqs pident) and `0.35` (ECFP Tanimoto),
respectively. With the default protein boundary, a normal below-boundary result such
as 12.5% is plotted as `0.125`; an at/above-boundary backfill is retained with a
consistency warning rather than being mislabeled as an expected low-similarity
recovery.

## Legacy CSV route

The old interface remains available for compatibility:

```text
--protein_training_data_path protein_training_data.csv
--ligand_training_data_path ligand_training_data.csv
```

These options provide already-derived CSVs. The new options provide reusable source
database bundles. Do not mix the two routes in one request. Legacy
`--build_bias_training_data` also remains confined to the legacy route.

If preparation reports a missing MMseqs executable, use
`cofolder-tools install-mmseqs` or set `COFOLDER_MMSEQS_BIN`. Incomplete downloads,
missing release metadata, stale checksums, wrong bundle kinds, and unsupported schema
versions are rejected with the affected file and a preparation command.
