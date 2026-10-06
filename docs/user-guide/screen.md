# Screen Command

For bias assessment during screening, see
[Providing bias training data](bias-training-data.md). Screen rows automatically
reuse invariant protein searches through the shared query cache.

The `screen` command runs **`validate` once per valid source record** in a CSV,
SDF, or MOL library. Malformed structure records are retained as failures while
valid records continue.

## Multi-parameter CSV screens

A CSV row can define a complete system variant by mapping any number of columns to
existing fields in the system YAML:

```bash
cofolder screen \
  -s system_screen.yaml -o options.yaml \
  --library parameter_screen.csv --col_id experiment \
  --map 'protein_sequence=sequences.0.protein.sequence' \
  --map 'ligand_smiles=sequences.1.ligand.smiles'
```

Each `--map` has the form `COLUMN=YAML_PATH`. Paths use dot-separated mapping keys
and zero-based list indices. Escape a literal dot or backslash in a key with a
backslash. Every destination must already exist in the template, all mapped cells
must be non-blank, and one CSV row always produces one system; COFOLDER does not
generate a Cartesian product.

String destinations preserve the CSV text exactly. Values mapped onto numbers,
booleans, lists, dictionaries, or `null` destinations are parsed as JSON and must
match the template type (`null` accepts any JSON value). This permits structured
values such as a complete constraints list to be supplied in one quoted CSV cell.

Mapped mode requires a CSV and `--col_id`, and it does not accept
`--smiles_column`. The latter remains part of the existing single-ligand library
mode. `--ligand_chain` is optional in mapped mode and selects the ligand for
ligand-specific filtering or clustering when those features are requested.

The completed system is validated independently for every row. A bad cell or system
fails only that row. If the template supplies an MSA and a mapping changes that
protein's sequence, the entire screen is rejected before prediction unless the MSA
is mapped with the sequence. This prevents a fixed alignment from being reused for
the wrong query:

```bash
cofolder screen \
  -s system_screen.yaml -o options.yaml \
  -c variants.csv --col_id experiment \
  --map 'protein_sequence=sequences.0.protein.sequence' \
  --map 'protein_msa=sequences.0.protein.msa'
```

Every mapped MSA is checked against its completed row sequence. Relative mapped MSA
paths resolve from the CSV directory. Alternatively, remove `msa` from the template;
COFOLDER will generate one alignment for each unique sequence and reuse it whenever
that sequence reappears. `msa: empty` remains an explicit request to skip generation.

Use `--preflight_only` to validate every completed row without running inference.
The report includes the mappings, valid and invalid row counts, row diagnostics, and
the total planned execution count.

## Basic Usage

```bash
cofolder screen \
  -s system.yaml \
  -o options.yaml \
  --library compounds.csv \
  --col_id compound_id \
  --ligand_chain B --smiles_column smiles
```

## Required Arguments

- `-s, --system_path`: Path to system YAML file
- `-o, --options_path`: Path to runner options YAML file
- `-c, --library`: Path to a CSV, SDF, or MOL library
- `--ligand_chain`: Existing ligand chain whose chemistry is replaced in
  single-ligand mode; optional analysis target in mapped mode

Single-ligand CSV libraries additionally require `--col_id` and
`--smiles_column`. Mapped CSV screens require `--col_id` and one or more `--map`
arguments. SDF/MOL libraries use `--id_property` (default `_Name`).

## Validation Rules

1. When supplied, `--ligand_chain` must resolve to exactly one ligand entity in
   every row.
2. Every source record receives a stable one-based `record_NNNNNN` identity.
3. Empty/invalid SMILES and malformed/unsanitizable molblocks fail only their record.
4. In single-ligand mode, protein, nucleic-acid, constraint, metadata, and
   unrelated-ligand fields remain fixed.

## Optional Arguments

- `-w, --wrk_dir`: Working directory
- `--merge_data`: Comma-separated metadata columns copied into summary
- `--library_format {csv,sdf,mol}`: Override extension-based format inference
- `--duplicate_id_policy {reject,suffix,source_index}`: Resolve duplicate IDs
  (`reject` by default)
- `--id_property`: SDF/MOL identifier property (`_Name` by default)
- `--map COLUMN=YAML_PATH`: Replace an existing system field from a CSV column;
  repeat to vary multiple parameters together
- All common validate options are supported and forwarded, including:
  - scoring functions
  - bias options (`--assess_bias`, `--bias_*`)
  - reproduction options (`--reference_path`, `--reproduction_metrics`)
  - robustness options
- `--ifp_filter_threshold FLOAT`: annotate compounds using inclusive reference
  overlap (`0` to `1`). Omit this option to disable filtering.
- `--ifp_filter_source {auto,reference_complex,custom_pocket}`: choose reference
  evidence; `auto` preserves custom-pocket precedence and otherwise uses
  `--reference_path`.
- `--ifp_taxonomy {distance,prolif}`: use normalized distance contacts or
  explicitly opt into subprocess-isolated ProLIF contacts.
- `--ifp_similarity_metric {jaccard,reference_coverage}`: choose the threshold
  score. Reference complexes default to Jaccard and custom pockets to coverage.
- `--ifp_filter_policy {similarity,required}` with repeatable
  `--ifp_required_interaction CHAIN:RESNUM[ICODE]:TYPE`: select threshold or
  required-contact behavior.
- `--ifp_reference_ligand CHAIN[:RESNUM[ICODE]]` and repeatable
  `--ifp_reference_receptor_chain CHAIN`: disambiguate reference entities.
- `--ligand_chain CHAIN` also selects the ligand evaluated by IFP filtering and clustering.
- `--cluster_ifps`: cluster evaluable typed IFPs after all rows finish.
- `--ifp_cluster_similarity_threshold FLOAT`: inclusive Jaccard-similarity cut
  for clustering (`0` to `1`, default `0.5`).

Custom-pocket filtering requires `--pocket_coverage_reference`; fingerprint
extraction is requested automatically. The reference can be a binary bitstring,
residue numbers or labels, or a file containing one of those forms:

```bash
cofolder screen \
  -s system.yaml -o options.yaml -c compounds.csv \
  --col_id compound_id \
  --ligand_chain B --smiles_column smiles \
  --scoring_functions ifp_distance \
  --pocket_coverage_reference "A25 G48 Y51" \
  --ifp_filter_threshold 0.6
```

Overlap is the fraction of active reference-pocket bits also present in the
predicted distance IFP. A row passes when overlap is greater than or equal to
the threshold. Filtering only annotates results; it never deletes rejected
rows or their prediction artifacts.

To derive the filter directly from a reference complex:

```bash
cofolder screen \
  -s system.yaml -o options.yaml -c compounds.csv \
  --col_id compound_id --ligand_chain B --smiles_column smiles \
  --reference_path reference_complex.cif \
  --ifp_filter_source reference_complex \
  --ifp_reference_ligand L:401 \
  --ifp_filter_threshold 0.6 \
  --ifp_similarity_metric jaccard
```

Unmappable structures remain in consolidated output with `not_evaluable`
status and machine-readable mapping diagnostics.

## IFP Clustering Without a Reference

Clustering does not require a reference structure or pocket definition. It is useful
for discovering recurring predicted contact patterns, including a reproducible
all-zero “no contacts” pattern:

```bash
cofolder screen \
  -s system.yaml -o options.yaml -c compounds.csv \
  --col_id compound_id \
  --ligand_chain B --smiles_column smiles \
  --scoring_functions ifp_distance \
  --cluster_ifps \
  --ifp_cluster_similarity_threshold 0.5
```

COFOLDER uses average-linkage agglomerative clustering over Jaccard distance and
cuts the tree at `distance <= 1 - similarity_threshold`. Stable IDs (`IFP001`,
`IFP002`, …) are assigned by each cluster's earliest input row. Missing, malformed,
or vector-length-incompatible IFPs receive `ifp_cluster_status=not_evaluable` and
no cluster ID. The first valid fingerprint defines the expected vector length.

## Full Diagnostic Screen

This command requests the structural, affinity, bias, distance-IFP, and pocket
diagnostics used for publication-facing review. Replace the reference and training
data paths with datasets appropriate for the target:

```bash
cofolder screen \
  -s system.yaml -o options.yaml -c compounds.csv \
  --col_id compound_id \
  --ligand_chain B --smiles_column smiles \
  --repeats 3 \
  --scoring_functions confidence_metrics affinity_metrics affinity_metrics_ext \
    sasa sasa_normalized ifp_distance \
  --assess_bias \
  --protein_training_data_path protein_training_data.csv \
  --ligand_training_data_path ligand_training_data.csv \
  --reference_path reference_complex.cif \
  --reproduction_metrics protein_rmsd ligand_rmsd sucos pocket_coverage \
  --pocket_coverage_reference "A25 G48 Y51" \
  --ifp_filter_threshold 0.6 \
  --cluster_ifps
```

If you only need pre-cofolding bias diagnostics for one system, prefer the dedicated
[`bias`](bias.md) workflow instead of running `screen`.

## Protein MSA Reuse

For Boltz-family runners, `screen` resolves each unique missing protein sequence once
and stores the reusable artifact under `<wrk_dir>/shared/msa/<runner>/`. The cache is
keyed by normalized sequence and MSA-generation settings, so a sequence order such as
`A, B, A, B` makes two MSA requests. The resolved MSA is injected into later repeats
and every subsequent row containing that sequence.

If the original system YAML already supplies `msa` for a protein, that file is used
directly and the server is not called for that protein. Relative MSA paths are resolved
relative to the original system YAML before row-specific YAML files are written. A
system may mix supplied and missing MSAs; only missing protein MSAs are generated.

In mapped mode, changing a sequence while inheriting such a supplied MSA is a setup
error. Map a matching MSA in the same row, or remove the template MSA to use automatic
sequence-keyed generation. Preflight reports the affected row, chain, mapping, and
alignment path before any MSA search or prediction starts.

The shared cache is matched by runner, protein sequence, and MSA-generation settings
(server URL, pairing strategy, maximum MSA depth, and backend version), and includes
a manifest. Credentials are not persisted. A missing, corrupt, sequence-mismatched,
or settings-incompatible cached artifact is not silently injected.

## Output

`screen` writes:

- Per-record validate outputs under `<wrk_dir>/compound_NNNNNN/...`
- Source/execution provenance under `<wrk_dir>/results/compound_members.csv`
- Mapping definitions under `<wrk_dir>/results/system_mappings.json` in mapped mode
- Canonical records: `<wrk_dir>/results/records.jsonl`
- Long-form views: `<wrk_dir>/results/{successes,metrics,failures}.csv`
- Manifest: `<wrk_dir>/results/manifest.json`
- Cluster summary, native linkage matrix, and leaf order (when enabled):
  `<wrk_dir>/results/ifp_cluster_{summary,linkage,leaf_order}.csv`
- Consolidated atom-level ProLIF occurrences (when ProLIF filtering or clustering is
  enabled): `<wrk_dir>/results/ifp_interaction_events.jsonl`

Summary columns include:

- `index`
- `<col_id>`
- `status` (`success` / `failed`)
- `error_message`
- `run_dir`
- the configured `smiles_column` value
- optional `merge_data` values
- filter audit columns: `ifp_filter_pass`, `ifp_filter_status`,
  `ifp_filter_reason`, `ifp_filter_overlap`, `ifp_filter_threshold`, and
  `ifp_filter_reference`, plus similarity, taxonomy, required-contact,
  missing-interaction, and mapping diagnostics
- clustering columns: `ifp_cluster_id` and `ifp_cluster_status`

Metric records include:

- all original input CSV columns
- run metadata (`index`, `status`, `error_message`, `run_dir`)
- extracted system-level score columns prefixed as `system__...`
- extracted chain-level score columns prefixed as `<entity>_<chain_id>__...`
- the same six filter audit columns as the summary CSV
- the same two clustering columns as the summary CSV

Filter status is `accepted`, `rejected`, `not_evaluable`, or `not_applied`.
Missing or malformed IFPs, all-zero references, and unequal vector lengths are
`not_evaluable`; vectors are never truncated for filtering.

The consolidated schema always includes the manuscript-facing confidence, affinity,
SASA, distance-IFP, training-set-proximity, and pocket-coverage columns. Unsupported or
unconfigured metrics are represented by empty values rather than omitted columns.
Pairwise confidence retains the runner's numeric-chain column and also exposes a
chain-ID alias such as `ligand_B__pair_chains_iptm_A`.

### Stable manuscript-facing score columns

The merged output always declares the following columns for configured system chain
IDs, using empty values when a backend or run does not provide a metric:

| Scope | Stable columns |
| --- | --- |
| System | `system__confidence_score`, `system__ptm`, `system__iptm` |
| Protein proximity | `system__bias_prot_sim_train_max`, `system__bias_prot_sim_train_pairwise_max`, `protein_<ID>__bias_prot_sim_train`, `protein_<ID>__bias_prot_sim_train_pairwise` |
| Ligand proximity | `system__bias_lig_sim_train_max`, `ligand_<ID>__bias_lig_sim_train` |
| Affinity/binding | `ligand_<ID>__affinity_pred_value`, `ligand_<ID>__affinity_probability_binary`, `ligand_<ID>__pIC50`, `ligand_<ID>__IC50_M`, `ligand_<ID>__pIC50_kcal_per_mol` |
| Chain confidence | `protein_<ID>__chains_ptm`, `ligand_<ID>__chains_ptm`, numeric `pair_chains_iptm_<index>` fields when emitted, and semantic `ligand_<ID>__pair_chains_iptm_<protein_ID>` aliases |
| Exposure and contacts | `ligand_<ID>__sasa`, `ligand_<ID>__sasa_norm_heavy`, `ligand_<ID>__ifp_distance`, `ligand_<ID>__ifp_prolif` |
| Pocket/reference | `system__pocket_coverage_ref`, `system__pocket_coverage_ref_mean`, `system__pocket_coverage_custom`, `system__pocket_coverage_custom_mean`, `ligand_<ID>__pocket_coverage_ref`, `ligand_<ID>__pocket_coverage_custom`, `ligand_<ID>__ligand_pose_overlap_ref` |

`ifp_cluster_status` is `clustered`, `not_evaluable`, or `not_applied`.
`ifp_cluster_summary.csv` contains `ifp_cluster_id`, `size`, JSON `member_ids`,
`medoid_compound_id`, JSON `consensus_ifp`, and
`mean_within_cluster_jaccard_similarity`. Consensus bits are present in at least
half of cluster members; medoid ties resolve to the earliest input row. A singleton
cluster has mean within-cluster similarity `1.0`.

`ifp_cluster_linkage.csv` is the native SciPy average-linkage matrix with explicit
merge indices. `ifp_cluster_leaf_order.csv` maps every dendrogram leaf back to its
input index, stable Screen execution key, and cluster ID. Empty and singleton runs
still publish these files with a complete header.

ProLIF fingerprints retain two representations. The compact `InteractionKey` set is
deduplicated by receptor residue and interaction type and remains the input to
filtering and clustering. The occurrence-level event records additionally retain
protein and ligand atom identities, ligand/protein roles, and named geometry.
Distances are in Å and angles are in degrees. PDB inputs are read directly; mmCIF
inputs are converted inside the isolated worker with original chain and residue
identities restored in the returned records.

The recorded linkage can be plotted without reclustering:

```python
from pathlib import Path

import pandas as pd
from scipy.cluster.hierarchy import dendrogram

results = Path("screen_output/results")
linkage = pd.read_csv(results / "ifp_cluster_linkage.csv")
leaves = pd.read_csv(results / "ifp_cluster_leaf_order.csv")
labels = leaves.sort_values("input_index")["member_id"].tolist()
dendrogram(
    linkage[["left_child", "right_child", "jaccard_distance", "member_count"]]
    .to_numpy(float),
    labels=labels,
)
```

Atom-specific classifications can be made from public event records. For example,
an Asp168 backbone-N contact is selected without invoking ProLIF directly:

```python
events = pd.read_json(results / "ifp_interaction_events.jsonl", lines=True)
asp168_n = events[events["protein_atoms"].apply(
    lambda atoms: any(
        atom["residue_name"] == "ASP"
        and atom["residue_number"] == 168
        and atom["atom_name"] == "N"
        for atom in atoms
    )
)]
```

When used from Python, `Screen.run()` returns the same merged results as a
`pandas.DataFrame` after writing the CSV files.

## Failure Behavior

If one row fails, screening continues for remaining rows.
Failures are recorded in `results/records.jsonl` and `results/failures.csv`.

## Related

- [Bias Command](bias.md)
- [Configuration Guide](../getting-started/configuration.md)
- [Virtual Screening Tutorial](../tutorials/screening.md)
- [Validate Command](validate.md)
