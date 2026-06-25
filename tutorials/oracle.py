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

    from _marimo_helpers import EXAMPLES_DIR, code_block, format_command, make_workspace, read_text, run_command

    return EXAMPLES_DIR, code_block, format_command, make_workspace, mo, read_text, run_command


@app.cell
def _(EXAMPLES_DIR, make_workspace):
    workspace = make_workspace("cofolder-tutorial-oracle-")
    system_path = EXAMPLES_DIR / "system.yaml"
    options_path = EXAMPLES_DIR / "options.yaml"
    return options_path, system_path, workspace


@app.cell
def _(code_block, format_command, mo, options_path, read_text, system_path, workspace):
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
        "affinity_pred_value",
        "--aggregate",
        "first",
        "-w",
        workspace,
        "--runner",
        "boltz2",
    ]
    mo.md(
        f"""
        # Oracle Workflow

        The `oracle` workflow wraps `validate` for one ligand input and returns one metric.

        **Options file**
        {code_block(read_text(options_path), "yaml")}

        **Command preview**
        ```bash
        {format_command(oracle_command)}
        ```
        """
    )
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
