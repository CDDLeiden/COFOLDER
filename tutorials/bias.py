import marimo

__generated_with = "0.23.8"
app = marimo.App()


@app.cell
def _():
    import sys
    from pathlib import Path
    import yaml

    notebook_dir = Path(__file__).resolve().parent
    if str(notebook_dir) not in sys.path:
        sys.path.insert(0, str(notebook_dir))

    import marimo as mo
    import pandas as pd

    from _marimo_helpers import code_block, format_command, make_workspace, read_text, run_command

    return code_block, format_command, make_workspace, mo, pd, read_text, run_command, yaml


@app.cell
def _(make_workspace, yaml):
    workspace = make_workspace("cofolder-tutorial-bias-")
    system_path = workspace / "bias_system.yaml"
    system_path.write_text(
        yaml.safe_dump(
            {
                "sequences": [
                    {"protein": {"id": "A", "sequence": "MKRAAT"}},
                    {"ligand": {"id": "B", "smiles": "CCO"}},
                ]
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    reference_structure = workspace / "custom_reference.pdb"
    reference_structure.write_text("HEADER CUSTOM\n", encoding="utf-8")
    custom_protein = workspace / "custom_protein.csv"
    custom_protein.write_text(
        "sequence,dataset_name,source_structure_path\n"
        f"MKRAAT,private_proteins,{reference_structure.name}\n",
        encoding="utf-8",
    )
    custom_ligand = workspace / "custom_ligand.csv"
    custom_ligand.write_text(
        "smiles,dataset_name,source_reference_path\n"
        f"CCO,private_ligands,{reference_structure.name}\n",
        encoding="utf-8",
    )
    return custom_ligand, custom_protein, system_path, workspace


@app.cell
def _(code_block, custom_ligand, custom_protein, format_command, mo, read_text, system_path, workspace):
    bias_command = [
        "cofolder",
        "bias",
        "--system_path",
        system_path,
        "--wrk_dir",
        workspace / "bias_run",
        "--custom_protein_reference_path",
        custom_protein,
        "--custom_ligand_reference_path",
        custom_ligand,
    ]
    mo.md(
        f"""
        # Bias Workflow

        This notebook walks through the standalone `bias` workflow using synthetic custom references.

        `cofolder bias` produces reference-overlap diagnostics for pre-cofolding decision support.
        These outputs are not binding-affinity predictions, model-confidence scores, or guarantees of
        cofolding success.

        **System file**
        {code_block(read_text(system_path), "yaml")}

        **Custom protein references**
        {code_block(read_text(custom_protein), "csv")}

        **Custom ligand references**
        {code_block(read_text(custom_ligand), "csv")}

        **Command preview**
        ```bash
        {format_command(bias_command)}
        ```
        """
    )
    return (bias_command,)


@app.cell
def _(mo):
    run_bias = mo.ui.checkbox(
        value=False,
        label="Run the standalone bias workflow in a temporary workspace",
    )
    run_bias
    return (run_bias,)


@app.cell
def _(bias_command, mo, pd, run_bias, run_command, workspace):
    mo.stop(
        not run_bias.value,
        mo.md("Enable the checkbox above to execute the bias tutorial."),
    )
    bias_output = run_command(bias_command, cwd=workspace)
    summary = pd.read_csv(
        workspace / "bias_run" / "results" / "bias_train" / "reference_landscape_summary.csv"
    )
    return bias_output, summary


@app.cell
def _(bias_output, mo, summary):
    mo.md(f"## Bias output\n```text\n{bias_output}\n```")
    summary
    return


if __name__ == "__main__":
    app.run()
