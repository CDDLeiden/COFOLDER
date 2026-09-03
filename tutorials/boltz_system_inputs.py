"""Interactive tutorial for DNA, RNA, and constraints with the Boltz family."""

import marimo

__generated_with = "0.23.8"
app = marimo.App()


@app.cell
def _():
    import tempfile
    from pathlib import Path

    import marimo as mo
    import pandas as pd

    from cofolder.acceptance.shared import (
        assert_backend_not_started,
        assert_chain_ids,
        assert_constraints_preserved,
        build_validate_command,
        combined_output,
        format_command,
        materialize_acceptance_inputs,
        materialize_invalid_constraint_input,
        options_resource_for_backend,
        reset_work_dir,
        run_cli_in_workspace,
    )
    from cofolder.modules.runners import get_runner

    return (
        Path,
        assert_backend_not_started,
        assert_chain_ids,
        assert_constraints_preserved,
        build_validate_command,
        combined_output,
        format_command,
        get_runner,
        materialize_acceptance_inputs,
        materialize_invalid_constraint_input,
        mo,
        options_resource_for_backend,
        pd,
        reset_work_dir,
        run_cli_in_workspace,
        tempfile,
    )


@app.cell
def _(mo):
    runner = mo.ui.dropdown(
        options=["boltz1", "boltz2", "boltz-community"],
        value="boltz2",
        label="Boltz runner",
    )
    run_prediction = mo.ui.checkbox(label="Run the real constrained prediction")
    run_rejection = mo.ui.checkbox(label="Run the preflight rejection example")
    mo.vstack([runner, run_prediction, run_rejection])
    return run_prediction, run_rejection, runner


@app.cell
def _(
    Path,
    build_validate_command,
    format_command,
    get_runner,
    materialize_acceptance_inputs,
    materialize_invalid_constraint_input,
    options_resource_for_backend,
    runner,
    tempfile,
):
    runner_name = runner.value
    workspace = Path(tempfile.mkdtemp(prefix=f"cofolder-tutorial-{runner_name}-inputs-"))
    fixtures = materialize_acceptance_inputs(
        workspace / "fixtures",
        options_resource=options_resource_for_backend(runner_name),
    )
    system_path = fixtures.constrained_system_paths[runner_name]
    runner_impl = get_runner(runner_name)
    runner_available, availability_message = runner_impl.check_availability()
    capabilities = runner_impl.input_capabilities

    prediction_dir = workspace / "prediction"
    prediction_command = build_validate_command(
        runner=runner_name,
        wrk_dir=prediction_dir,
        system_path=system_path,
        options_path=fixtures.options_path,
        scoring_functions=["confidence_metrics"],
    )
    invalid_path = materialize_invalid_constraint_input(
        runner_name,
        system_path,
        workspace / "fixtures" / f"invalid_{runner_name}.yaml",
    )
    rejection_dir = workspace / "rejection"
    rejection_command = build_validate_command(
        runner=runner_name,
        wrk_dir=rejection_dir,
        system_path=invalid_path,
        options_path=fixtures.options_path,
        scoring_functions=["confidence_metrics"],
    )
    install_command = f'python -m pip install -e ".[tutorials,{runner_name}]"'
    return (
        availability_message,
        capabilities,
        install_command,
        prediction_command,
        prediction_dir,
        rejection_command,
        rejection_dir,
        runner_available,
        runner_name,
        system_path,
        workspace,
    )


@app.cell
def _(
    availability_message,
    capabilities,
    format_command,
    install_command,
    mo,
    prediction_command,
    runner_available,
    runner_name,
    system_path,
    workspace,
):
    source_yaml = system_path.read_text(encoding="utf-8")
    availability = availability_message or (
        f"Runner '{runner_name}' is available."
        if runner_available
        else f"Runner '{runner_name}' is not installed."
    )
    mo.md(
        f"""
        # DNA, RNA, and constraints with the Boltz family

        This tutorial sends one small **protein → DNA → RNA → ligand** system through
        COFOLDER's real `validate` workflow. Residue numbers in contacts are **1-based**.
        Boltz-1 accepts bonds and one 6 Å pocket; Boltz-2 and boltz-community also accept
        contact constraints.

        Install the selected backend in its own environment:

        ```bash
        {install_command}
        marimo edit tutorials/boltz_system_inputs.py
        ```

        **Backend status:** `{availability}`

        **Declared entity types:** `{', '.join(sorted(capabilities.entity_types))}`<br>
        **Declared constraint types:** `{', '.join(sorted(capabilities.constraint_types))}`

        The real run requires a CUDA-capable machine and may download backend data. It is
        disabled until you select the checkbox above. Outputs stay in `{workspace}`.

        ## Canonical COFOLDER input

        ```yaml
        {source_yaml.rstrip()}
        ```

        ## Command preview

        ```bash
        {format_command(prediction_command)}
        ```
        """
    )
    return


@app.cell
def _(
    assert_chain_ids,
    assert_constraints_preserved,
    mo,
    pd,
    prediction_command,
    prediction_dir,
    reset_work_dir,
    run_cli_in_workspace,
    run_prediction,
    runner_available,
    system_path,
    workspace,
):
    mo.stop(not run_prediction.value, "Enable the real constrained prediction checkbox.")
    mo.stop(
        not runner_available,
        "Install the selected Boltz backend in this environment before running inference.",
    )
    reset_work_dir(prediction_dir)
    prediction_result = run_cli_in_workspace(prediction_command, workspace)
    prepared_path = prediction_dir / "raw" / system_path.name
    preserved_constraints = assert_constraints_preserved(system_path, prepared_path)
    chain_metrics_path = prediction_dir / "results" / "chain_metrics.csv"
    assert_chain_ids(chain_metrics_path, ["A", "D", "R", "L"])
    chain_metrics = pd.read_csv(chain_metrics_path)
    system_metrics = pd.read_csv(prediction_dir / "results" / "system_metrics.csv")
    prepared_yaml = prepared_path.read_text(encoding="utf-8")
    return (
        chain_metrics,
        prediction_result,
        prepared_yaml,
        preserved_constraints,
        system_metrics,
    )


@app.cell
def _(chain_metrics, mo, prepared_yaml, preserved_constraints, system_metrics):
    mo.vstack(
        [
            mo.md(
                f"""
                ## Successful backend run

                COFOLDER preserved all `{len(preserved_constraints)}` top-level constraint(s)
                while preparing the ligand. The YAML below is the actual system file submitted
                to Boltz.

                ```yaml
                {prepared_yaml.rstrip()}
                ```

                The ordered metadata proves that protein, DNA, RNA, and ligand chains all
                survived normalization.
                """
            ),
            mo.md("### Chain metrics"),
            mo.ui.table(chain_metrics),
            mo.md("### System metrics"),
            mo.ui.table(system_metrics),
        ]
    )
    return


@app.cell
def _(
    assert_backend_not_started,
    combined_output,
    mo,
    rejection_command,
    rejection_dir,
    reset_work_dir,
    run_cli_in_workspace,
    run_rejection,
    runner_available,
    runner_name,
    workspace,
):
    mo.stop(not run_rejection.value, "Enable the preflight rejection checkbox.")
    mo.stop(
        not runner_available,
        "Install the selected Boltz backend before running the CLI rejection example.",
    )
    reset_work_dir(rejection_dir)
    rejection_result = run_cli_in_workspace(
        rejection_command,
        workspace,
        check=False,
    )
    if rejection_result.returncode == 0:
        raise AssertionError("The intentionally invalid input was unexpectedly accepted.")
    rejection_output = combined_output(rejection_result)
    expected_error = (
        "type 'contact' is unsupported"
        if runner_name == "boltz1"
        else "references unknown chain 'Z'"
    )
    if expected_error not in rejection_output:
        raise AssertionError(
            f"Expected validation message {expected_error!r}, received:\n{rejection_output}"
        )
    assert_backend_not_started(rejection_dir)
    mo.md(
        f"""
        ## Rejected before inference

        Exit code: `{rejection_result.returncode}`. No `raw/repeat_*` directory was created,
        so the backend was never launched.

        ```text
        {rejection_output[-4000:]}
        ```
        """
    )
    return


if __name__ == "__main__":
    app.run()
