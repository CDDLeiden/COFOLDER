# Metric reference

This table mirrors `cofolder.modules.contracts.METRIC_CATALOG`, the authoritative
registry used to validate public records. A metric appears only when its scope and
evidence are applicable. Direction describes optimization semantics, not a universal
scientific threshold.

| Metric | Group | Class | Type | Unit/scale | Direction | Scope | Evidence |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `confidence_score` | `confidence_metrics` | confidence | float | unitless | maximize | chain, system | reference_free |
| `ptm` | `confidence_metrics` | confidence | float | unitless | maximize | chain, system | reference_free |
| `iptm` | `confidence_metrics` | confidence | float | unitless | maximize | chain, system | reference_free |
| `chains_ptm` | `confidence_metrics` | confidence | float | unitless | maximize | chain, system | reference_free |
| `chain_ptm` | `confidence_metrics` | confidence | float | unitless | maximize | chain, system | reference_free |
| `sample_ranking_score` | `confidence_metrics` | confidence | float | unitless | maximize | chain, system | reference_free |
| `avg_plddt` | `confidence_metrics` | confidence | float | unitless | maximize | chain, system | reference_free |
| `ligand_iptm` | `confidence_metrics` | confidence | float | unitless | maximize | system | reference_free |
| `protein_iptm` | `confidence_metrics` | confidence | float | unitless | maximize | system | reference_free |
| `complex_plddt` | `confidence_metrics` | confidence | float | unitless | maximize | system | reference_free |
| `complex_iplddt` | `confidence_metrics` | confidence | float | unitless | maximize | system | reference_free |
| `complex_pde` | `confidence_metrics` | confidence | float | unitless | minimize | system | reference_free |
| `complex_ipde` | `confidence_metrics` | confidence | float | unitless | minimize | system | reference_free |
| `complex_pae` | `confidence_metrics` | confidence | float | unitless | minimize | system | reference_free |
| `complex_ipae` | `confidence_metrics` | confidence | float | unitless | minimize | system | reference_free |
| `gpde` | `confidence_metrics` | confidence | float | unitless | minimize | chain, system | reference_free |
| `disorder` | `confidence_metrics` | confidence | float | unitless | minimize | chain, system | reference_free |
| `has_clash` | `confidence_metrics` | confidence | boolean | unitless | minimize | chain, system | reference_free |
| `pair_chains_iptm` | `confidence_metrics` | confidence | float | unitless | maximize | chain_pair | reference_free |
| `chain_pair_iptm` | `confidence_metrics` | confidence | float | unitless | maximize | chain_pair | reference_free |
| `bespoke_iptm` | `confidence_metrics` | confidence | float | unitless | maximize | chain_pair | reference_free |
| `affinity_pred_value` | `affinity_metrics` | affinity | float | log10(µM) | minimize | chain | reference_free |
| `affinity_probability_binary` | `affinity_metrics` | binding | float | unitless | maximize | chain | reference_free |
| `pIC50` | `affinity_metrics_ext` | affinity | float | pIC50 | maximize | chain | reference_free |
| `IC50_M` | `affinity_metrics_ext` | affinity | float | M | minimize | chain | reference_free |
| `IC50_uM` | `affinity_metrics_ext` | affinity | float | µM | minimize | chain | reference_free |
| `pIC50_kcal_per_mol` | `affinity_metrics_ext` | affinity | float | kcal/mol | maximize | chain | reference_free |
| `sasa` | `structure_metrics` | structural | float | Å² | neutral | chain | reference_free |
| `sasa_norm_heavy` | `structure_metrics` | structural | float | Å²/heavy atom | neutral | chain | reference_free |
| `ifp_distance` | `structure_metrics` | structural | JSON | unitless | neutral | chain | reference_free |
| `ifp_distance_features` | `structure_metrics` | structural | JSON | unitless | neutral | chain | reference_free |
| `ifp_prolif` | `structure_metrics` | structural | JSON | unitless | neutral | chain | reference_free |
| `ifp_prolif_features` | `structure_metrics` | structural | JSON | unitless | neutral | chain | reference_free |
| `ligand_rmsd_ref` | `reproduction_metrics` | reproduction | float | Å | minimize | chain, system | reference_structure |
| `protein_rmsd_ref` | `reproduction_metrics` | reproduction | float | Å | minimize | chain, system | reference_structure |
| `ligand_rmsd_ref_mean` | `reproduction_metrics` | reproduction | float | Å | minimize | chain, system | reference_structure |
| `protein_rmsd_ref_mean` | `reproduction_metrics` | reproduction | float | Å | minimize | chain, system | reference_structure |
| `sucos_ref` | `reproduction_metrics` | reproduction | float | unitless | maximize | chain, system | reference_structure |
| `sucos_shape_ref` | `reproduction_metrics` | reproduction | float | unitless | maximize | chain, system | reference_structure |
| `sucos_feature_ref` | `reproduction_metrics` | reproduction | float | unitless | maximize | chain, system | reference_structure |
| `sucos_ref_mean` | `reproduction_metrics` | reproduction | float | unitless | maximize | chain, system | reference_structure |
| `ligand_pose_overlap_ref` | `reproduction_metrics` | reproduction | float | unitless | maximize | chain, system | reference_structure |
| `pocket_coverage_ref` | `reproduction_metrics` | reproduction | float | unitless | maximize | chain, system | reference_structure |
| `pocket_coverage_ref_mean` | `reproduction_metrics` | reproduction | float | unitless | maximize | chain, system | reference_structure |
| `pocket_coverage_custom` | `reproduction_metrics` | reproduction | float | unitless | maximize | chain, system | custom_pocket |
| `pocket_coverage_custom_mean` | `reproduction_metrics` | reproduction | float | unitless | maximize | chain, system | custom_pocket |
| `bias_prot_sim_train` | `bias_metrics` | training_proximity | float | percent identity | minimize | chain, system | reference_free |
| `bias_prot_sim_train_max` | `bias_metrics` | training_proximity | float | percent identity | minimize | chain, system | reference_free |
| `bias_prot_sim_train_pairwise` | `bias_metrics` | training_proximity | float | percent identity | minimize | chain, system | reference_free |
| `bias_prot_sim_train_pairwise_max` | `bias_metrics` | training_proximity | float | percent identity | minimize | chain, system | reference_free |
| `bias_lig_sim_train` | `bias_metrics` | training_proximity | float | Tanimoto | minimize | chain, system | reference_free |
| `bias_lig_sim_train_max` | `bias_metrics` | training_proximity | float | Tanimoto | minimize | chain, system | reference_free |
| `ifp_filter_pass` | `screen_metrics` | filter | boolean | unitless | neutral | compound | custom_pocket, reference_structure |
| `ifp_filter_overlap` | `screen_metrics` | filter | float | unitless | maximize | compound | custom_pocket, reference_structure |
| `ifp_filter_similarity` | `screen_metrics` | filter | float | unitless | maximize | compound | custom_pocket, reference_structure |
| `ifp_filter_threshold` | `screen_metrics` | filter | float | unitless | neutral | compound | custom_pocket, reference_structure |
| `ifp_filter_status` | `screen_metrics` | filter | string | unitless | neutral | compound | custom_pocket, reference_structure |
| `ifp_filter_reason` | `screen_metrics` | filter | string | unitless | neutral | compound | custom_pocket, reference_structure |
| `ifp_filter_reference` | `screen_metrics` | filter | string | unitless | neutral | compound | custom_pocket, reference_structure |
| `ifp_filter_similarity_metric` | `screen_metrics` | filter | string | unitless | neutral | compound | custom_pocket, reference_structure |
| `ifp_filter_policy` | `screen_metrics` | filter | string | unitless | neutral | compound | custom_pocket, reference_structure |
| `ifp_filter_taxonomy` | `screen_metrics` | filter | string | unitless | neutral | compound | custom_pocket, reference_structure |
| `ifp_filter_mapping_status` | `screen_metrics` | filter | string | unitless | neutral | compound | custom_pocket, reference_structure |
| `ifp_filter_required_interactions` | `screen_metrics` | filter | JSON | unitless | neutral | compound | custom_pocket, reference_structure |
| `ifp_filter_missing_interactions` | `screen_metrics` | filter | JSON | unitless | neutral | compound | custom_pocket, reference_structure |
| `ifp_filter_mapping_failures` | `screen_metrics` | filter | JSON | unitless | neutral | compound | custom_pocket, reference_structure |
| `ifp_cluster_id` | `screen_metrics` | filter | string | unitless | neutral | compound | reference_free |
| `ifp_cluster_status` | `screen_metrics` | filter | string | unitless | neutral | compound | reference_free |
| `oracle_score` | `oracle_metrics` | oracle | float | unitless | neutral | compound | all regimes |
| `oracle_raw_score` | `oracle_metrics` | oracle | float | unitless | neutral | compound | all regimes |
| `oracle_gate_adjusted_score` | `oracle_metrics` | oracle | float | unitless | neutral | compound | all regimes |
| `struct_rmsd` | `robustness_metrics` | robustness | float | Å | minimize | chain, system | reference_free |

## Availability and missing values

Every metric record has one of `computed`, `missing`, `unsupported`, `failed`, or
`not_requested`. Only `computed` records contain a finite value; every other state
uses JSON `null` plus a machine-readable reason when applicable. Confidence,
binding likelihood, model-derived affinity, structural agreement, robustness, and
training proximity remain distinct diagnostic classes. A direction is not a claim
that one metric or threshold is universally valid.
