import marimo

__generated_with = "0.23.8"
app = marimo.App()


@app.cell
def _():
    import sys
    from pathlib import Path

    notebook_dir = Path(__file__).resolve().parent
    if str(notebook_dir) not in sys.path:
        sys.path.insert(0, str(notebook_dir))

    import marimo as mo

    from _marimo_helpers import (
        EXAMPLES_DIR,
        REPO_ROOT,
        code_block,
        format_command,
        make_workspace,
        read_text,
        run_command,
    )

    return (
        EXAMPLES_DIR,
        REPO_ROOT,
        code_block,
        format_command,
        make_workspace,
        mo,
        read_text,
        run_command,
    )


@app.cell
def _(EXAMPLES_DIR, REPO_ROOT, make_workspace):
    workspace = make_workspace("cofolder-tutorial-oracle-")
    system_path = EXAMPLES_DIR / "system.yaml"
    options_path = EXAMPLES_DIR / "options.yaml"
    ligand_training_path = (
        REPO_ROOT / "src/cofolder/acceptance/data/ligand_training_data.csv"
    )
    return ligand_training_path, options_path, system_path, workspace


@app.cell
def _(
    code_block,
    format_command,
    ligand_training_path,
    mo,
    options_path,
    read_text,
    system_path,
    workspace,
):
    oracle_command = [
        "cofolder",
        "oracle",
        "-s",
        system_path,
        "-o",
        options_path,
        "--input_smiles",
        "CCO",
        "--output_metric",
        "ligand_B__affinity_pred_value",
        "--aggregate",
        "first",
        "-w",
        workspace,
        "--runner",
        "boltz2",
    ]
    bias_command = [
        "cofolder",
        "oracle",
        "-s",
        system_path,
        "-o",
        options_path,
        "--input_smiles",
        "CCO",
        "--output_metric",
        "ligand_B__bias_lig_sim_train",
        "--assess_bias",
        "--bias_chains",
        "B",
        "--ligand_training_data_path",
        ligand_training_path,
        "-w",
        workspace / "bias",
    ]
    pocket_command = [
        "cofolder",
        "oracle",
        "-s",
        system_path,
        "-o",
        options_path,
        "--input_smiles",
        "CCO",
        "--scoring_functions",
        "ifp_distance",
        "--reproduction_metrics",
        "pocket_coverage",
        "--pocket_coverage_reference",
        "A25 G48 Y51",
        "--output_metric",
        "ligand_B__pocket_coverage_custom",
        "-w",
        workspace / "pocket",
    ]
    python_example = "\n".join(
        [
            "from cofolder.recipes.oracle import Oracle, OracleGate, OracleGatePolicy",
            "",
            "oracle = Oracle(",
            '    wrk_dir="oracle_gated",',
            '    system_path="examples/system.yaml",',
            '    options_path="examples/options.yaml",',
            '    input_smiles="CCO",',
            '    aggregate="mean",',
            "    scoring_functions=[",
            '        "confidence_metrics", "affinity_metrics_ext",',
            '        "ifp_distance", "sasa_normalized",',
            "    ],",
            '    pocket_coverage_reference="A25 G48 Y51",',
            '    reproduction_metrics=["pocket_coverage"],',
            "    score_components={",
            '        "ligand_B__pIC50": 1.0,',
            '        "system__confidence_score": 1.0,',
            "    },",
            "    score_gates=[",
            '        OracleGate("ligand_B__pocket_coverage_custom", "ge", 0.60),',
            '        OracleGate("ligand_B__sasa_norm_heavy", "le", 2.0),',
            "    ],",
            '    gate_policy=OracleGatePolicy("downweight", 0.25),',
            ")",
            "score = oracle.run()",
        ]
    )
    mo.md(f"""
        # Oracle Workflow

        The `oracle` workflow wraps `validate` for one ligand input and returns one
        scalar. Qualified names make the metric's system or chain scope explicit.

        **Options file**
        {code_block(read_text(options_path), "yaml")}

        ## Affinity command
        ```bash
        {format_command(oracle_command)}
        ```

        ## Ligand training-set similarity
        ```bash
        {format_command(bias_command)}
        ```

        ## Custom pocket/IFP-reference coverage
        ```bash
        {format_command(pocket_command)}
        ```

        ## Python composite with structural gates
        {code_block(python_example, "python")}

        The thresholds above are illustrative. Calibrate pocket and SASA gates during
        validation for the target system. Replace the down-weight policy with
        `OracleGatePolicy("fixed_penalty", value)` or
        `OracleGatePolicy("non_binder", value)` when that interpretation fits the
        external optimization workflow.
        """)
    return (oracle_command,)


@app.cell
def _(mo):
    run_oracle = mo.ui.checkbox(
        value=False,
        label="Run the oracle workflow in a temporary workspace",
    )
    run_oracle
    return (run_oracle,)


@app.cell
def _(mo, oracle_command, run_command, run_oracle, workspace):
    mo.stop(
        not run_oracle.value,
        mo.md("Enable the checkbox above to execute the oracle tutorial."),
    )
    oracle_output = run_command(oracle_command, cwd=workspace)
    return (oracle_output,)


@app.cell
def _(mo, oracle_output):
    mo.md(f"## Command output\n```text\n{oracle_output}\n```")
    return


if __name__ == "__main__":
    app.run()
