"""Interactive tutorial for DNA, RNA, and pocket constraints with OpenFold3."""

import marimo

__generated_with = "0.23.8"
app = marimo.App()


@app.cell
def _():
    import json
    import os
    import tempfile
    from pathlib import Path

    import marimo as mo
    import pandas as pd

    from cofolder.acceptance.shared import (
        assert_backend_not_started,
        assert_chain_ids,
        assert_constraints_preserved,
        assert_openfold3_query_translation,
        build_validate_command,
        combined_output,
        format_command,
        inspect_openfold3_notebook_setup,
        materialize_acceptance_inputs,
        materialize_invalid_constraint_input,
        options_resource_for_backend,
        reset_work_dir,
        resolve_openfold3_notebook_cache,
        rewrite_openfold3_options_cache_path,
        run_cli_in_workspace,
    )
    from cofolder.modules.runners import get_runner

    return (
        Path,
        assert_backend_not_started,
        assert_chain_ids,
        assert_constraints_preserved,
        assert_openfold3_query_translation,
        build_validate_command,
        combined_output,
        format_command,
        get_runner,
        inspect_openfold3_notebook_setup,
        json,
        materialize_acceptance_inputs,
        materialize_invalid_constraint_input,
        mo,
        options_resource_for_backend,
        os,
        pd,
        reset_work_dir,
        resolve_openfold3_notebook_cache,
        rewrite_openfold3_options_cache_path,
        run_cli_in_workspace,
        tempfile,
    )


@app.cell
def _(mo, os):
    cache_path = mo.ui.text(
        value=os.environ.get("OPENFOLD_CACHE", ""),
        label="OpenFold3 cache path",
        full_width=True,
    )
    run_prediction = mo.ui.checkbox(label="Run the real constrained prediction")
    run_rejection = mo.ui.checkbox(label="Run the unsupported-bond rejection")
    mo.vstack([cache_path, run_prediction, run_rejection])
    return cache_path, run_prediction, run_rejection


@app.cell
def _(
    Path,
    build_validate_command,
    cache_path,
    format_command,
    get_runner,
    inspect_openfold3_notebook_setup,
    materialize_acceptance_inputs,
    materialize_invalid_constraint_input,
    options_resource_for_backend,
    resolve_openfold3_notebook_cache,
    rewrite_openfold3_options_cache_path,
    tempfile,
):
    runner_name = "openfold3"
    workspace = Path(tempfile.mkdtemp(prefix="cofolder-tutorial-openfold3-inputs-"))
    fixtures = materialize_acceptance_inputs(
        workspace / "fixtures",
        options_resource=options_resource_for_backend(runner_name),
    )
    system_path = fixtures.constrained_system_paths[runner_name]
    runner_impl = get_runner(runner_name)
    runner_available, availability_message = runner_impl.check_availability()
    capabilities = runner_impl.input_capabilities

    cache_config = resolve_openfold3_notebook_cache(override=cache_path.value)
    setup_status = inspect_openfold3_notebook_setup(cache_config)
    if cache_config.configured_cache is not None:
        rewrite_openfold3_options_cache_path(
            fixtures.options_path,
            cache_config.configured_cache,
        )

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
        workspace / "fixtures" / "invalid_openfold3_bond.yaml",
    )
    rejection_dir = workspace / "rejection"
    rejection_command = build_validate_command(
        runner=runner_name,
        wrk_dir=rejection_dir,
        system_path=invalid_path,
        options_path=fixtures.options_path,
        scoring_functions=["confidence_metrics"],
    )
    return (
        availability_message,
        cache_config,
        capabilities,
        prediction_command,
        prediction_dir,
        rejection_command,
        rejection_dir,
        runner_available,
        runner_name,
        setup_status,
        system_path,
        workspace,
    )


@app.cell
def _(
    availability_message,
    capabilities,
    format_command,
    mo,
    prediction_command,
    runner_available,
    setup_status,
    system_path,
    workspace,
):
    source_yaml = system_path.read_text(encoding="utf-8")
    availability = availability_message or (
        "OpenFold3 is available." if runner_available else "OpenFold3 is not installed."
    )
    mo.md(
        f"""
        # DNA, RNA, and pocket constraints with OpenFold3

        COFOLDER accepts one canonical system schema across runners. This tutorial shows how
        its protein, DNA, RNA, ligand, and pocket data become a real OpenFold3 query. Pocket
        residue numbers are **1-based**. OpenFold3 accepts one pocket constraint; bond and
        contact constraints are rejected because this adapter cannot safely pass them through.

        Prepare a dedicated environment and cache before opening the notebook:

        ```bash
        python -m pip install -e ".[tutorials,openfold3]"
        export OPENFOLD_CACHE="$PWD/cache/.openfold3-cache"
        scripts/setup_openfold3.sh
        marimo edit tutorials/openfold3_system_inputs.py
        ```

        **Backend status:** `{availability}`<br>
        **Setup status:** `{setup_status.message}`<br>
        **Declared entity types:** `{', '.join(sorted(capabilities.entity_types))}`<br>
        **Declared constraint types:** `{', '.join(sorted(capabilities.constraint_types))}`

        The real run requires a CUDA-capable machine and is disabled until selected. Outputs
        stay in `{workspace}`.

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
    assert_openfold3_query_translation,
    cache_config,
    json,
    mo,
    pd,
    prediction_command,
    prediction_dir,
    reset_work_dir,
    run_cli_in_workspace,
    run_prediction,
    runner_available,
    setup_status,
    system_path,
    workspace,
):
    mo.stop(not run_prediction.value, "Enable the real constrained prediction checkbox.")
    mo.stop(not runner_available, "Install OpenFold3 before running inference.")
    mo.stop(not setup_status.ready, setup_status.message)
    reset_work_dir(prediction_dir)
    prediction_result = run_cli_in_workspace(
        prediction_command,
        workspace,
        env=cache_config.env,
    )
    prepared_path = prediction_dir / "raw" / system_path.name
    preserved_constraints = assert_constraints_preserved(system_path, prepared_path)
    query_path = (
        prediction_dir
        / "raw"
        / "repeat_1"
        / "openfold3"
        / "inference_query_set.json"
    )
    translated_query = assert_openfold3_query_translation(
        system_path,
        query_path,
        system_name=system_path.stem,
    )
    chain_metrics_path = prediction_dir / "results" / "chain_metrics.csv"
    assert_chain_ids(chain_metrics_path, ["A", "D", "R", "L"])
    chain_metrics = pd.read_csv(chain_metrics_path)
    system_metrics = pd.read_csv(prediction_dir / "results" / "system_metrics.csv")
    query_json = json.dumps({"queries": {system_path.stem: translated_query}}, indent=2)
    return (
        chain_metrics,
        prediction_result,
        preserved_constraints,
        query_json,
        system_metrics,
    )


@app.cell
def _(chain_metrics, mo, preserved_constraints, query_json, system_metrics):
    mo.vstack(
        [
            mo.md(
                f"""
                ## Successful backend run

                COFOLDER preserved `{len(preserved_constraints)}` canonical pocket constraint
                and translated it to OpenFold3's query-level `pocket_constraint`. This is the
                actual query JSON consumed by OpenFold3:

                ```json
                {query_json}
                ```
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
    cache_config,
    combined_output,
    mo,
    rejection_command,
    rejection_dir,
    reset_work_dir,
    run_cli_in_workspace,
    run_rejection,
    runner_available,
    workspace,
):
    mo.stop(not run_rejection.value, "Enable the unsupported-bond rejection checkbox.")
    mo.stop(not runner_available, "Install OpenFold3 before running the CLI rejection example.")
    reset_work_dir(rejection_dir)
    rejection_result = run_cli_in_workspace(
        rejection_command,
        workspace,
        check=False,
        env=cache_config.env,
    )
    if rejection_result.returncode == 0:
        raise AssertionError("The intentionally unsupported bond was unexpectedly accepted.")
    rejection_output = combined_output(rejection_result)
    expected_error = "type 'bond' is unsupported"
    if expected_error not in rejection_output:
        raise AssertionError(
            f"Expected validation message {expected_error!r}, received:\n{rejection_output}"
        )
    assert_backend_not_started(rejection_dir)
    mo.md(
        f"""
        ## Rejected before inference

        Exit code: `{rejection_result.returncode}`. No `raw/repeat_*` directory was created,
        proving that COFOLDER rejected the unsupported bond before launching OpenFold3.

        ```text
        {rejection_output[-4000:]}
        ```
        """
    )
    return


if __name__ == "__main__":
    app.run()
