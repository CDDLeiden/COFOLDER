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
    workspace = make_workspace("cofolder-tutorial-screen-")
    system_path = EXAMPLES_DIR / "system_screen.yaml"
    options_path = EXAMPLES_DIR / "options.yaml"
    screen_csv = EXAMPLES_DIR / "ligand_screen.csv"
    return options_path, screen_csv, system_path, workspace


@app.cell
def _(code_block, format_command, mo, options_path, read_text, screen_csv, system_path, workspace):
    screen_command = [
        "cofolder",
        "screen",
        "-s",
        system_path,
        "-o",
        options_path,
        "-c",
        screen_csv,
        "--col_id",
        "Name",
        "--variable",
        "sequences,1,ligand,smiles",
        "--col_variable",
        "SMILES",
        "--merge_data",
        "pIC50",
        "-w",
        workspace,
        "--runner",
        "boltz2",
        "--cluster_ifps",
    ]
    mo.md(
        f"""
        # Screen Workflow

        This notebook covers the `screen` workflow, which applies `validate` row by row.

        **System template**
        {code_block(read_text(system_path), "yaml")}

        **Options file**
        {code_block(read_text(options_path), "yaml")}

        **Screening CSV**
        {code_block(read_text(screen_csv), "csv")}

        **Command preview**
        ```bash
        {format_command(screen_command)}
        ```
        """
    )
    return (screen_command,)


@app.cell
def _(mo):
    run_screen = mo.ui.checkbox(
        value=False,
        label="Run the screen workflow in a temporary workspace",
    )
    run_screen
    return (run_screen,)


@app.cell
def _(mo, run_command, run_screen, screen_command, workspace):
    mo.stop(
        not run_screen.value,
        mo.md("Enable the checkbox above to execute the screening tutorial."),
    )
    screen_output = run_command(screen_command, cwd=workspace)
    return (screen_output,)


@app.cell
def _(mo, screen_output):
    mo.md(f"## Command output\n```text\n{screen_output}\n```")
    return


if __name__ == "__main__":
    app.run()
