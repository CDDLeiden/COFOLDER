import marimo

__generated_with = "0.23.8"
app = marimo.App(width="medium")


@app.cell
def _():
    import sys
    from pathlib import Path

    notebook_dir = Path(__file__).resolve().parent
    if str(notebook_dir) not in sys.path:
        sys.path.insert(0, str(notebook_dir))

    import marimo as mo
    import pandas as pd

    from _marimo_helpers import code_block, make_workspace, read_text
    from _structure_gated_oracle import (
        ASSET_DIR,
        CONFIG_PATH,
        ILLUSTRATIVE_RESULTS_PATH,
        analyze_pose,
        fetch_references,
        load_config,
        load_ligands,
        make_oracle_scoring_function,
        molecule_grid_svg,
        rank_candidates,
        reanalyze_saved_predictions,
        reference_geometry_report,
        write_audit_bundle,
    )

    return (
        ASSET_DIR,
        CONFIG_PATH,
        ILLUSTRATIVE_RESULTS_PATH,
        analyze_pose,
        code_block,
        fetch_references,
        load_config,
        load_ligands,
        make_oracle_scoring_function,
        make_workspace,
        mo,
        molecule_grid_svg,
        pd,
        rank_candidates,
        reanalyze_saved_predictions,
        read_text,
        reference_geometry_report,
        write_audit_bundle,
    )


@app.cell
def _(ASSET_DIR, load_config, load_ligands, make_workspace, molecule_grid_svg):
    policy = load_config()
    ligands = load_ligands()
    workspace = make_workspace("cofolder-mapk14-structure-gated-")
    ligand_svg = molecule_grid_svg(ligands)
    return ligand_svg, ligands, policy, workspace


@app.cell
def _(ligand_svg, ligands, mo):
    panel_intro = mo.md(
        f"""
        # MAPK14 structure-gated Oracle

        This tutorial implements a strict ranking rule:

        **prefer Type II binding → prefer more designated interactions → use the
        predicted affinity score to break ties.**

        Literature annotations are expectations used to design this validation
        panel. They never enter the computed reward. The pose predicted for each
        candidate supplies every structural decision.

        Install the notebook, analysis, and Boltz2 dependencies with:

        ```bash
        python -m pip install -e ".[analysis,tutorials,boltz2]"
        marimo edit tutorials/structure_gated_oracle.py
        ```

        ## Six-ligand panel

        {ligand_svg}
        """
    )
    panel_table = mo.ui.table(
        ligands[
            [
                "candidate_id",
                "name",
                "role",
                "literature_endpoint",
                "literature_value_nm",
                "literature_source",
            ]
        ],
        selection=None,
    )
    mo.vstack([panel_intro, panel_table])
    return


@app.cell
def _(CONFIG_PATH, code_block, mo, read_text):
    mo.md(
        f"""
        ## Frozen structural policy

        Type II requires contact with both the ATP site and back pocket plus the
        manuscript SI D1/D2 DFG-out rule. Three literature-grounded hydrogen-bond
        features are counted once
        each. The threshold is fixed from 1KV2 and 1A9U before candidate predictions
        are inspected.

        {code_block(read_text(CONFIG_PATH), "yaml")}

        The scalar is `(K + 1) * is_type_II + matched_count + bounded(raw_score)`,
        where the bounded term is the arctangent CDF.
        Boltz2 reports `affinity_pred_value` as lower-is-better log10(µM), so the
        callback negates it before applying this higher-is-better rule.
        Because the transformed score is strictly between zero and one, affinity
        cannot overturn one interaction and interactions cannot overturn Type II
        eligibility.
        """
    )
    return


@app.cell
def _(ILLUSTRATIVE_RESULTS_PATH, pd, rank_candidates):
    illustrative_input = pd.read_csv(ILLUSTRATIVE_RESULTS_PATH)
    illustrative_ranked = rank_candidates(illustrative_input)
    return illustrative_input, illustrative_ranked


@app.cell
def _(illustrative_ranked, mo):
    slide_table = illustrative_ranked[
        [
            "candidate_id",
            "structural_classification",
            "dfg_d1_angstrom",
            "dfg_d2_angstrom",
            "glu71_sidechain_hbond",
            "asp168_backbone_hbond",
            "met109_backbone_hbond",
            "matched_interactions",
            "native_affinity_score",
            "raw_score",
            "scalar_reward",
            "score_rank",
            "structure_gated_rank",
        ]
    ].sort_values("structure_gated_rank")
    walkthrough_intro = mo.md(
        """
        ## Inexpensive walkthrough

        These rows are explicitly illustrative. They exercise the same ranking code
        as a live run, but they are not Boltz2 predictions and are not experimental
        measurements. Notice that BIRB796 rises above the higher raw score of
        compound 48 because it matches one more designated interaction.
        """
    )
    slide_widget = mo.ui.table(slide_table, selection=None)
    mo.vstack([walkthrough_intro, slide_widget])
    return (slide_table,)


@app.cell
def _(
    CONFIG_PATH,
    ILLUSTRATIVE_RESULTS_PATH,
    illustrative_input,
    mo,
    workspace,
    write_audit_bundle,
):
    illustrative_csv, illustrative_evidence, illustrative_provenance = (
        write_audit_bundle(
            illustrative_input,
            workspace / "illustrative_audit",
            inputs=(CONFIG_PATH, ILLUSTRATIVE_RESULTS_PATH),
            metadata={"data_kind": "illustrative", "seed": None},
        )
    )
    mo.md(
        f"""
        The complete illustrative audit was written to
        `{illustrative_csv}`, atom-level evidence in `{illustrative_evidence}`, and
        provenance in `{illustrative_provenance}`.
        """
    )
    return illustrative_csv, illustrative_evidence, illustrative_provenance


@app.cell
def _(mo):
    check_references = mo.ui.checkbox(
        value=False,
        label="Download 1KV2 and 1A9U and verify the manuscript DFG D1/D2 rule",
    )
    check_references
    return (check_references,)


@app.cell
def _(check_references, fetch_references, mo, reference_geometry_report, workspace):
    mo.stop(not check_references.value, mo.md("Reference download is optional."))
    reference_paths = fetch_references(workspace / "references")
    geometry_report = reference_geometry_report(reference_paths)
    geometry_report
    return geometry_report, reference_paths


@app.cell
def _(mo, policy):
    key_rows = "\n".join(
        f"- `{feature['key']}` on receptor atom(s) "
        f"`{', '.join(feature['receptor_atoms'])}` — {feature['rationale']}"
        for feature in policy["interaction_policy"]["features"]
    )
    mo.md(
        f"""
        ## What is inspected in each predicted complex

        The pose must contact at least one configured ATP-site residue and at least
        one back-pocket residue. The manuscript SI rule additionally requires
        Asn155 CA–Phe169 CA D1 ≤
        `{policy["structural_predicate"]["dfg_conformation"]["d1"]["threshold_angstrom"]}` Å
        and Glu71 CA–Phe169 CA D2 ≥
        `{policy["structural_predicate"]["dfg_conformation"]["d2"]["threshold_angstrom"]}` Å.
        After curated ligand graph recovery and PDB2PQR protein hydrogen placement,
        ProLIF tests these ligand-centric features using its frozen 3.5 Å and
        130–180° hydrogen-bond defaults:

        {key_rows}

        Open a predicted mmCIF in PyMOL or ChimeraX and highlight MAPK14 residues
        `51, 53, 71, 75, 84, 86, 104, 106, 109, 168, 169`. The downloaded reference
        structures provide a direct Type I/Type II comparison.
        """
    )
    return


@app.cell
def _(mo):
    saved_run_root = mo.ui.text(
        value="",
        label="Saved six-ligand run directory",
        placeholder="/path/to/completed/run",
    )
    reanalyse_saved = mo.ui.checkbox(
        value=False,
        label="Reanalyse saved structures and exact backend affinity records",
    )
    mo.vstack(
        [
            mo.md(
                "## Reanalyse saved predictions\n\nThis path performs no GPU work. "
                "Each score is read from its run's `metrics.csv` and kept paired "
                "with that run's single saved structure."
            ),
            saved_run_root,
            reanalyse_saved,
        ]
    )
    return reanalyse_saved, saved_run_root


@app.cell
def _(
    mo,
    policy,
    rank_candidates,
    reanalyse_saved,
    reanalyze_saved_predictions,
    saved_run_root,
):
    mo.stop(
        not reanalyse_saved.value or not saved_run_root.value.strip(),
        mo.md("Enter a completed run directory and enable reanalysis."),
    )
    saved_input = reanalyze_saved_predictions(
        saved_run_root.value.strip(), config=policy
    )
    saved_ranked = rank_candidates(saved_input, config=policy)
    mo.ui.table(
        saved_ranked.sort_values("structure_gated_rank", na_position="last"),
        selection=None,
    )
    return saved_input, saved_ranked


@app.cell
def _(mo):
    run_live = mo.ui.checkbox(
        value=False,
        label="Run all six MAPK14 ligands with Boltz2 (GPU and service access required)",
    )
    run_live
    return (run_live,)


@app.cell
def _(
    ASSET_DIR,
    analyze_pose,
    ligands,
    make_oracle_scoring_function,
    mo,
    pd,
    policy,
    run_live,
    workspace,
):
    mo.stop(
        not run_live.value,
        mo.md("Enable the checkbox only when the Boltz2 backend is ready."),
    )
    from cofolder.recipes.oracle import Oracle

    live_rows = []
    live_failures = []
    for ligand in ligands.itertuples(index=False):
        callback = make_oracle_scoring_function(
            candidate_id=ligand.candidate_id,
            expected_role=ligand.role,
            canonical_smiles=ligand.canonical_smiles,
            audit_sink=live_rows,
            config=policy,
        )
        try:
            Oracle(
                wrk_dir=str(workspace / "live" / ligand.candidate_id),
                system_path=str(ASSET_DIR / "system.yaml"),
                options_path=str(ASSET_DIR / "options.yaml"),
                runner="boltz2",
                input_smiles=ligand.canonical_smiles,
                ligand_chain="B",
                repeats=1,
                seed=int(policy["execution"]["seed"]),
                scoring_functions=["affinity_metrics"],
                assess_robustness=False,
                scoring_function=callback,
            ).run()
        except Exception as exc:
            if not any(row["candidate_id"] == ligand.candidate_id for row in live_rows):
                live_rows.append(
                    analyze_pose(
                        workspace / "live" / ligand.candidate_id / "missing.cif",
                        raw_score=float("nan"),
                        candidate_id=ligand.candidate_id,
                        prediction_id="unavailable",
                        expected_role=ligand.role,
                        config=policy,
                    )
                )
            live_failures.append(
                {"candidate_id": ligand.candidate_id, "error": str(exc)}
            )
    live_input = pd.DataFrame(live_rows)
    return live_failures, live_input


@app.cell
def _(
    CONFIG_PATH,
    ASSET_DIR,
    live_failures,
    live_input,
    mo,
    rank_candidates,
    workspace,
    write_audit_bundle,
):
    live_ranked = rank_candidates(live_input)
    live_csv, live_evidence, live_provenance = write_audit_bundle(
        live_input,
        workspace / "live_audit",
        inputs=(
            CONFIG_PATH,
            ASSET_DIR / "ligands.csv",
            ASSET_DIR / "system.yaml",
            ASSET_DIR / "options.yaml",
        ),
        metadata={
            "data_kind": "boltz2_prediction",
            "seed": 20260928,
            "failures": live_failures,
        },
    )
    live_intro = mo.md(
        f"""
        ## Live result

        Audit: `{live_csv}`  
        Atom-level evidence: `{live_evidence}`  
        Provenance: `{live_provenance}`

        Literature roles remain visible for comparison but have no effect on either
        reward or rank. A failed structure, interaction extraction, or score is
        retained as `not_evaluable` and receives no rank.
        """
    )
    live_table = mo.ui.table(
        live_ranked.sort_values("structure_gated_rank", na_position="last"),
        selection=None,
    )
    mo.vstack([live_intro, live_table])
    return live_csv, live_evidence, live_provenance, live_ranked


@app.cell
def _(code_block, mo):
    reuse = """audit_rows = []
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
        input_smiles="YOUR_SMILES",
        ligand_chain="B",
        repeats=1,
        scoring_functions=["affinity_metrics"],
        assess_robustness=False,
        scoring_function=callback,
    ).run()
    """
    mo.md(
        f"""
        ## Reuse the scalar Oracle

        {code_block(reuse, "python")}

        Keep one repeat and one diffusion sample unless you first define a policy
        that keeps pose, interaction evidence, and score paired within every sample.
        """
    )
    return


if __name__ == "__main__":
    app.run()
