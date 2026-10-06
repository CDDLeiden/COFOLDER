# MAPK14 Structure-gated Oracle

This tutorial turns a predicted MAPK14–ligand complex into a scalar while preserving
a strict scientific priority:

1. prefer a Type II pose;
2. among poses in the same tier, prefer more designated interactions;
3. use the model affinity score only to break the remaining ties.

The executable version is `tutorials/structure_gated_oracle.py`. It contains an
offline walkthrough and an optional six-ligand Boltz2 run. Install and open it with:

```bash
python -m pip install -e ".[analysis,tutorials,boltz2]"
marimo edit tutorials/structure_gated_oracle.py
```

## What the example does

The fixed target is canonical human p38α/MAPK14 ([UniProt Q16539](https://www.uniprot.org/uniprotkb/Q16539)).
The same unrestrained protein input and prediction settings are used for every
candidate. Reference structures define the analysis policy; they are not supplied to
Boltz2 as pocket restraints.

| Candidate | Literature role | Evidence used to construct the panel |
| --- | --- | --- |
| SB203580 | Type I inhibitor | ATP-site pyridinyl-imidazole; Type I contrast structure [PDB 1A9U](https://www.rcsb.org/structure/1A9U) |
| Compound 43 | lower-affinity Type II | BIRB796 matched series, calculated Kd 14 nM |
| Compound 48 | intermediate-affinity Type II | same series and assay, calculated Kd 0.52 nM |
| BIRB796, compound 5 | higher-affinity Type II | same series and assay, calculated Kd 0.046 nM; [PDB 1KV2](https://www.rcsb.org/structure/1KV2) |
| SB202474 | expected inactive control | established inactive p38 comparator; inactivity is not proof of no binding |
| UM101 | expected Type IV/remote-site binder | experimental binding with a computationally proposed substrate-docking site; no co-crystal pose |

Compounds 43, 48, and 5 come from one thermal-denaturation study and use the same
calculated-Kd endpoint. Compound 48 is the tested member closest to the logarithmic
midpoint between 43 and 5. The study describes their shared DFG-out scaffold and the
ethoxy-linked substituents directed toward the ATP region. See
[Regan et al. 2003](https://doi.org/10.1021/jm030121k) and the preceding
[discovery study](https://doi.org/10.1021/jm020057r).

The SB202474 label comes from its use as an inactive comparator in
[primary literature](https://pmc.ncbi.nlm.nih.gov/articles/PMC2171902/). UM101's
binding and proposed ED-site model come from
[Shah et al.](https://doi.org/10.4049/jimmunol.1602059). These labels provide
expectations for validation. They are never inputs to the computed reward.

Canonical structures, endpoints, citations, and evidence notes are versioned in
`tutorials/assets/structure_gated_oracle/ligands.csv`.

## Freeze the structural rule first

The policy is in `tutorials/assets/structure_gated_oracle/config.yaml`. A predicted
pose qualifies as Type II only when all three conditions hold:

- a ligand heavy atom lies within 4.5 Å of an ATP-site residue;
- a ligand heavy atom lies within 4.5 Å of a back-pocket residue;
- the COFOLDER manuscript SI classical DFG-out test passes: D1, the Asn155
  CA–Phe169 CA distance, is at most 7.2 Å, and D2, the Glu71 CA–Phe169 CA
  distance, is at least 9.0 Å.

The Type II 1KV2 reference gives D1 = 5.949 Å and D2 = 11.076 Å. The Type I 1A9U
reference gives D1 = 8.096 Å and D2 = 8.555 Å. The notebook can download the two
mmCIF files, verify their pinned SHA-256 hashes, and recompute both distances.
Residue numbers refer to the canonical UniProt sequence and predicted protein chain
`A`.

## Prepare chemistry before measuring interactions

The pose file supplies coordinates, while each curated isomeric SMILES supplies the
ligand molecular graph: connectivity, bond orders, aromaticity, formal charges, and
specified stereochemistry. The tutorial maps that graph onto the pose's heavy atoms,
adds ligand hydrogens without moving heavy atoms, and records the atom correspondence.
An identity or mapping disagreement makes the candidate `not_evaluable`.

The protein is prepared with PDB2PQR's AMBER templates at pH 7.4. Hydrogen placement
is optimized, histidine assignments and any added missing heavy atoms are recorded,
and every supplied protein heavy-atom coordinate is restored and checked. ProLIF runs
in a separate worker process on these prepared molecules. The baseline preserves the
protonation state encoded by each curated ligand SMILES; protonation states are not
selected from the resulting ranks.

The three features use **ligand-centric ProLIF terminology**:

```text
ligand donor   → Glu71 side-chain OE1/OE2   A:71:hb_donor
ligand acceptor ← Asp168 backbone N         A:168:hb_acceptor
ligand acceptor ← Met109 backbone N         A:109:hb_acceptor
```

The Glu71 and Asp168 urea interactions are described in the
[BIRB796 discovery paper](https://doi.org/10.1021/jm020057r). The
[matched-series SAR paper](https://doi.org/10.1021/jm030121k) motivates testing the
additional Met109 backbone interaction made by the ATP-region tail. ProLIF's installed
hydrogen-bond defaults are frozen at a 3.5 Å donor–acceptor distance and a 130–180°
D–H–A angle. A residue/type event counts only when its receptor atom and donor/acceptor
roles match the table. Duplicate atom-level events still count once. Hydrophobic,
aromatic, ionic, and van der Waals contacts remain in the JSON diagnostics and do not
contribute to reward.

The reference check recovers all three published interactions in 1KV2: Glu71 at
2.977 Å and 159.72°, Asp168 at 2.890 Å and 163.29°, and Met109 at 2.900 Å and
152.14°. The same preparation applied to 1A9U recovers its Met109 hinge interaction;
1A9U is a DFG-in structural contrast and is not assumed to lack hinge binding.

These definitions are a worked MAPK14 policy, not universal kinase criteria.
Recalibrate residue sets and geometry before applying the tutorial to another target.

## Turn the rule into one scalar

Let `K` be the number of designated interactions. The callback returns:

```text
(K + 1) × is_type_II + matched_count + bounded(raw_score)
```

The bounded term uses `0.5 + atan(raw_score) / π`. It is strictly increasing and
lies between zero and one.
Boltz2's `affinity_pred_value` is a lower-is-better log10(µM) quantity, so the
callback multiplies it by `-1`. The tutorial's `raw_score` is this oriented,
higher-is-better value.
It therefore preserves the affinity ordering without allowing affinity to overturn
one interaction. With `K + 1` as the Type II weight, all interactions combined
cannot overturn Type II eligibility.

This scalar has the same order as sorting descending by:

```text
(is_type_II, matched_count, raw_score)
```

The notebook's illustrative table exercises this exact code and shows both raw-score
and structure-gated ranks. Those rows are teaching data, not experimental results or
Boltz2 predictions.

## Run the live simulation

The live checkbox runs `Oracle` once for each ligand. Every call uses one repeat and
one diffusion sample so the structure, interaction fingerprint, and affinity score
come from the same prediction:

```python
from cofolder.recipes.oracle import Oracle
from _structure_gated_oracle import make_oracle_scoring_function

audit_rows = []
callback = make_oracle_scoring_function(
    candidate_id="my_candidate",
    expected_role="unknown prospectively",
    canonical_smiles="YOUR_SMILES",
    audit_sink=audit_rows,
)

reward = Oracle(
    wrk_dir="oracle_my_candidate",
    system_path="tutorials/assets/structure_gated_oracle/system.yaml",
    options_path="tutorials/assets/structure_gated_oracle/options.yaml",
    runner="boltz2",
    input_smiles="YOUR_SMILES",
    ligand_chain="B",
    repeats=1,
    seed=20260928,
    scoring_functions=["affinity_metrics"],
    assess_robustness=False,
    scoring_function=callback,
).run()
```

Do not increase repeats or diffusion samples without defining a pose-aware ensemble
policy. Independently aggregating affinity, geometry, and interactions could combine
evidence from different poses.

## Inspect and audit the result

Load each predicted mmCIF in PyMOL, ChimeraX, or another molecular viewer. Compare it
with 1A9U and 1KV2. This PyMOL setup gives consistent molecular views of the
configured sites and the two SI distances:

```pml
select atp_site, chain A and resi 51+53+106+109
select back_pocket, chain A and resi 71+75+84+86+104+168+169
select um101_region, chain A and resi 177+180+181+182+183+186+252+255+292+294
distance dfg_d1, chain A and resi 155 and name CA, chain A and resi 169 and name CA
distance dfg_d2, chain A and resi 71 and name CA, chain A and resi 169 and name CA
show sticks, organic or atp_site or back_pocket or um101_region
```

The notebook also renders the six ligand structures and reports the reference
geometry and interaction checks.
Ready-to-run views are provided as `view_interactions.pml` and
`view_interactions.cxc` in the tutorial asset directory.

To inspect an existing run, enter its root directory in the notebook's saved-run
section. The reanalysis reads the six structures and exact affinity records and does
not invoke Boltz2. The expected layout is
`<root>/<candidate_id>/oracle_run/results/{structures,metrics.csv}`.

The audit CSV records candidate and prediction identity, all three feature states,
D1/D2, Type II eligibility, the native affinity value, its oriented score, scalar
reward, both ranks, and evaluation status. The interaction JSON sidecar retains ligand
and receptor atom identities, distances, angles, all diagnostic contacts, preparation
settings, charges, atom mapping, heavy-coordinate checks, and package versions. The
provenance JSON records the seed, policy hash, and input and structure hashes.

A missing structure, failed interaction extraction, or unavailable score is recorded
as `not_evaluable` and receives no rank. It is not interpreted as a non-binder.
Predictions are allowed to disagree with the literature labels; that disagreement is
one of the validation outcomes this panel is designed to reveal.

The six saved Boltz2 predictions from 2026-09-28 were reanalysed without rerunning the
GPU backend. Exact affinity values came from each saved `metrics.csv`, rather than
rounded console output. All six were evaluable. The observed gated order was BIRB796,
compound 48, compound 43, UM101, SB203580, then SB202474. BIRB796 and compound 48
each matched all three features, so oriented affinity broke their tie. Compound 43
matched Glu71 and Asp168 but missed Met109. The predicted UM101 pose unexpectedly
matched Glu71 and Asp168 in the kinase pocket; it did not pass the Type II predicate,
and this disagreement with its literature expectation remains visible rather than
being forced away.

The machine-readable record is
`tutorials/assets/structure_gated_oracle/acceptance.json`, its summary is
`acceptance-v2.csv`, and atom-level evidence is in
`acceptance-v2-interactions.json`. The earlier four-feature policy is archived as
`acceptance-v1.json` and marked superseded. Full atom-level reference recovery and
preparation audits are retained in `reference-validation.json`.

## Related

- [Oracle user guide](../user-guide/oracle.md)
- [Virtual screening tutorial](screening.md)
- [Metric reference](../user-guide/metric-reference.md)
