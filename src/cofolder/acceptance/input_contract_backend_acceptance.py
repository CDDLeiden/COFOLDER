"""Opt-in Marimo lane for real-backend system-input acceptance."""

import marimo

__generated_with = "0.23.8"
app = marimo.App()


@app.cell
def _():
    import os
    import tempfile
    from pathlib import Path

    import marimo as mo

    from cofolder.acceptance.shared import (
        assert_chain_ids,
        assert_runner_available,
        assert_runner_setup_ready,
        build_validate_command,
        materialize_acceptance_inputs,
        options_resource_for_backend,
        reset_work_dir,
        rewrite_openfold3_options_cache_path,
        run_cli_in_workspace,
        scoring_functions_for_backend_acceptance,
    )

    return (
        Path,
        assert_chain_ids,
        assert_runner_available,
        assert_runner_setup_ready,
        build_validate_command,
        materialize_acceptance_inputs,
        mo,
        options_resource_for_backend,
        os,
        reset_work_dir,
        rewrite_openfold3_options_cache_path,
        run_cli_in_workspace,
        scoring_functions_for_backend_acceptance,
        tempfile,
    )


@app.cell
def _(mo):
    runner = mo.ui.dropdown(
        options=["boltz1", "boltz2", "boltz-community", "openfold3"],
        value="boltz2",
        label="Runner",
    )
    run_nucleic = mo.ui.checkbox(label="Run DNA/RNA acceptance")
    run_constraint = mo.ui.checkbox(label="Run constraint acceptance")
    mo.vstack([runner, run_nucleic, run_constraint])
    return run_constraint, run_nucleic, runner


@app.cell
def _(
    Path,
    assert_runner_available,
    assert_runner_setup_ready,
    materialize_acceptance_inputs,
    options_resource_for_backend,
    os,
    rewrite_openfold3_options_cache_path,
    runner,
    tempfile,
):
    runner_name = runner.value
    workspace = Path(tempfile.mkdtemp(prefix=f"cofolder-{runner_name}-input-contract-"))
    fixtures = materialize_acceptance_inputs(
        workspace / "fixtures",
        options_resource=options_resource_for_backend(runner_name),
    )
    assert_runner_available(runner_name)
    env = dict(os.environ)
    if runner_name == "openfold3":
        assert_runner_setup_ready(runner_name, env=env)
        cache_path = Path(env.get("OPENFOLD_CACHE", "~/.openfold3")).expanduser()
        rewrite_openfold3_options_cache_path(fixtures.options_path, cache_path)
    return env, fixtures, runner_name, workspace


@app.cell
def _(
    assert_chain_ids,
    build_validate_command,
    env,
    fixtures,
    mo,
    reset_work_dir,
    run_cli_in_workspace,
    run_nucleic,
    runner_name,
    scoring_functions_for_backend_acceptance,
    workspace,
):
    mo.stop(not run_nucleic.value, "Enable the DNA/RNA acceptance checkbox.")
    output = reset_work_dir(workspace / "nucleic")
    run_cli_in_workspace(
        build_validate_command(
            runner=runner_name,
            wrk_dir=output,
            system_path=fixtures.nucleic_acid_system_path,
            options_path=fixtures.options_path,
            scoring_functions=scoring_functions_for_backend_acceptance(runner_name),
        ),
        workspace,
        env=env,
    )
    assert_chain_ids(output / "results" / "chain_metrics.csv", ["A", "D", "R", "L"])
    return


@app.cell
def _(
    assert_chain_ids,
    build_validate_command,
    env,
    fixtures,
    mo,
    reset_work_dir,
    run_cli_in_workspace,
    run_constraint,
    runner_name,
    scoring_functions_for_backend_acceptance,
    workspace,
):
    mo.stop(not run_constraint.value, "Enable the constraint acceptance checkbox.")
    output = reset_work_dir(workspace / "constraint")
    run_cli_in_workspace(
        build_validate_command(
            runner=runner_name,
            wrk_dir=output,
            system_path=fixtures.constrained_system_paths[runner_name],
            options_path=fixtures.options_path,
            scoring_functions=scoring_functions_for_backend_acceptance(runner_name),
        ),
        workspace,
        env=env,
    )
    assert_chain_ids(output / "results" / "chain_metrics.csv", ["A", "D", "R", "L"])
    return


if __name__ == "__main__":
    app.run()
