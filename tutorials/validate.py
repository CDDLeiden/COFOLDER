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

    from _marimo_helpers import EXAMPLES_DIR, code_block, format_command, make_workspace, read_text, run_command
    from _workflow_examples import build_validate_command

    return EXAMPLES_DIR, build_validate_command, code_block, format_command, make_workspace, mo, pd, read_text, run_command


@app.cell
def _(EXAMPLES_DIR, build_validate_command, make_workspace):
    workspace = make_workspace("cofolder-tutorial-validate-")
    system_path = EXAMPLES_DIR / "system.yaml"
    options_path = EXAMPLES_DIR / "options.yaml"
    validate_command = build_validate_command(EXAMPLES_DIR, workspace)
    return options_path, system_path, validate_command, workspace


@app.cell
def _(code_block, format_command, mo, options_path, read_text, system_path, validate_command):
    mo.md(
        f"""
        # Validate Workflow

        This notebook walks through the single-system `validate` workflow.

        **System file**
        {code_block(read_text(system_path), "yaml")}

        **Options file**
        {code_block(read_text(options_path), "yaml")}

        **Command preview**
        ```bash
        {format_command(validate_command)}
        ```
        """
    )
    return (validate_command,)


@app.cell
def _(mo):
    run_validate = mo.ui.checkbox(
        value=False,
        label="Run the validate workflow in a temporary workspace",
    )
    run_validate
    return (run_validate,)


@app.cell
def _(mo, run_command, run_validate, validate_command, workspace):
    mo.stop(
        not run_validate.value,
        mo.md("Enable the checkbox above to execute the validate tutorial."),
    )
    validate_output = run_command(validate_command, cwd=workspace)
    return (validate_output,)


@app.cell
def _(mo, validate_output):
    mo.md(f"## Command output\n```text\n{validate_output}\n```")
    return


if __name__ == "__main__":
    app.run()
