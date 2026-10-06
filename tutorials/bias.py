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
    import pandas as pd

    from _marimo_helpers import code_block, format_command, make_workspace, read_text, run_command
    from _workflow_examples import build_bias_command, prepare_bias_example

    return (
        build_bias_command,
        code_block,
        format_command,
        make_workspace,
        mo,
        pd,
        prepare_bias_example,
        read_text,
        run_command,
    )


@app.cell
def _(build_bias_command, make_workspace, prepare_bias_example):
    workspace = make_workspace("cofolder-tutorial-bias-")
    inputs = prepare_bias_example(workspace)
    bias_command = build_bias_command(inputs, workspace / "bias_run")
    return bias_command, inputs, workspace


@app.cell
def _(bias_command, code_block, format_command, inputs, mo, read_text):
    mo.md(
        f"""
        # Bias Workflow

        This notebook walks through the standalone `bias` workflow using synthetic custom references.

        `cofolder bias` produces reference-overlap diagnostics for pre-cofolding decision support.
        These outputs are not binding-affinity predictions, model-confidence scores, or guarantees of
        cofolding success.

        **System file**
        {code_block(read_text(inputs.system_path), "yaml")}

        **Custom protein references**
        {code_block(read_text(inputs.protein_reference_path), "csv")}

        **Custom ligand references**
        {code_block(read_text(inputs.ligand_reference_path), "csv")}

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
