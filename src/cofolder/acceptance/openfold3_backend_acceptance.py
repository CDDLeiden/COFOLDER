import marimo

__generated_with = "0.23.8"
app = marimo.App()


@app.cell
def _():
    import os
    import subprocess
    import tempfile
    from pathlib import Path

    import marimo as mo

    from cofolder.acceptance.shared import (
        assert_runner_setup_ready,
        assert_csv_columns_have_values,
        assert_file_exists,
        assert_runner_available,
        build_oracle_command,
        build_screen_command,
        build_validate_command,
        format_command,
        install_command_for_backend,
        materialize_acceptance_inputs,
        options_resource_for_backend,
        oracle_metric_for_runner,
        oracle_scoring_functions_for_runner,
        reset_work_dir,
        scoring_functions_for_backend_acceptance,
    )

    def stream_cli_in_notebook(
        command: list[str],
        workspace: Path,
        *,
        env: dict[str, str] | None = None,
    ) -> str:
        print(f"$ {format_command(command)}", flush=True)
        process = subprocess.Popen(
            command,
            cwd=str(workspace),
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        if process.stdout is None:
            raise RuntimeError("Expected subprocess stdout pipe to be available for live logging.")

        lines: list[str] = []
        for line in process.stdout:
            print(line, end="", flush=True)
            lines.append(line)

        return_code = process.wait()
        output = "".join(lines)
        if return_code != 0:
            raise RuntimeError(
                f"Command failed with exit code {return_code}: {format_command(command)}\n"
                f"STREAMED OUTPUT:\n{output}"
            )
        return output

    return (
        Path,
        assert_csv_columns_have_values,
        assert_file_exists,
        assert_runner_available,
        assert_runner_setup_ready,
        build_oracle_command,
        build_screen_command,
        build_validate_command,
        format_command,
        install_command_for_backend,
        materialize_acceptance_inputs,
        mo,
        options_resource_for_backend,
        oracle_metric_for_runner,
        oracle_scoring_functions_for_runner,
        os,
        reset_work_dir,
        scoring_functions_for_backend_acceptance,
        stream_cli_in_notebook,
        tempfile,
    )


@app.cell
def _(mo):
    mo.md("""
    # OpenFold3 backend acceptance

    Recommended fresh-environment install path for this notebook:

    ```bash
    python -m pip install -e ".[acceptance,openfold3]"
    ```

    This is a dedicated OpenFold3 lane. Before running the expensive cells:

    - prepare the OpenFold3 cache and model data with `scripts/setup_openfold3.sh`
    - keep `OPENFOLD_CACHE` pointed at a deliberate location for this environment when you want something other than the standard `~/.openfold3`
    - treat this lane as confidence-only: it requests `confidence_metrics`, not affinity groups
    - keep MSA/template experiments separate unless you are deliberately testing them

    Live log behavior:

    - each command streams stdout/stderr into the cell console while it runs
    - the final captured output is also shown in the rendered notebook output after completion
    """)
    return


@app.cell
def _(mo):
    run_validate = mo.ui.checkbox(value=False, label="Run validate acceptance step")
    run_validate
    return (run_validate,)


@app.cell
def _(mo):
    run_screen = mo.ui.checkbox(value=False, label="Run screen acceptance step")
    run_screen
    return (run_screen,)


@app.cell
def _(mo):
    run_oracle = mo.ui.checkbox(value=False, label="Run oracle acceptance step")
    run_oracle
    return (run_oracle,)


@app.cell
def _(
    Path,
    assert_runner_available,
    install_command_for_backend,
    materialize_acceptance_inputs,
    options_resource_for_backend,
    os,
    tempfile,
):
    runner = "openfold3"
    workspace = Path(tempfile.mkdtemp(prefix="cofolder-openfold3-acceptance-"))
    fixtures = materialize_acceptance_inputs(
        workspace / "fixtures",
        options_resource=options_resource_for_backend(runner),
    )
    openfold_cache = workspace / ".openfold3"
    openfold_cache.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env["OPENFOLD_CACHE"] = str(openfold_cache)
    install_command = install_command_for_backend(runner)
    availability_message = assert_runner_available(runner)
    return (
        availability_message,
        env,
        fixtures,
        install_command,
        openfold_cache,
        runner,
        workspace,
    )


@app.cell
def _(availability_message, install_command, mo, openfold_cache):
    mo.md(f"""
    **Verified install command:** `{install_command}`

    **Availability check:** `{availability_message}`

    **Setup requirement:** run `scripts/setup_openfold3.sh` against the cache path below before starting any expensive workflow cell.

    **OPENFOLD_CACHE for this notebook:** `{openfold_cache}`

    Before triggering the expensive steps below, run upstream setup in the same environment, for example:

    ```bash
    export OPENFOLD_CACHE="{openfold_cache}"
    scripts/setup_openfold3.sh
    ```

    The OpenFold3 extra installs the backend package, while the setup script prepares the cache, checkpoints, and CCD in the chosen cache root.
    """)
    return


@app.cell
def _(
    fixtures,
    format_command,
    mo,
    runner,
    scoring_functions_for_backend_acceptance,
    workspace,
):
    scoring = scoring_functions_for_backend_acceptance(runner)
    mo.md(
        f"""
        **Workspace:** `{workspace}`

        **OPENFOLD_CACHE:** `{workspace / ".openfold3"}`

        **Validate scoring groups:** `{', '.join(scoring)}`

        **Validate command preview**
        ```bash
        {format_command([
            "cofolder", "validate", "-s", str(fixtures.system_path), "-o", str(fixtures.options_path),
            "-w", str(workspace / "validate"), "--runner", runner, "--repeats", "1", "--scoring_functions", *scoring
        ])}
        ```
        """
    )
    return (scoring,)


@app.cell
def _(
    assert_csv_columns_have_values,
    assert_file_exists,
    assert_runner_setup_ready,
    build_validate_command,
    env,
    fixtures,
    mo,
    reset_work_dir,
    run_validate,
    runner,
    scoring,
    stream_cli_in_notebook,
    workspace,
):
    mo.stop(
        not run_validate.value,
        mo.md("Enable `Run validate acceptance step` above before starting the expensive validate run."),
    )
    assert_runner_setup_ready(runner, env=env)
    validate_dir = reset_work_dir(workspace / "validate")
    validate_command = build_validate_command(
        runner=runner,
        wrk_dir=validate_dir,
        system_path=fixtures.system_path,
        options_path=fixtures.options_path,
        scoring_functions=scoring,
    )
    validate_output = stream_cli_in_notebook(validate_command, workspace, env=env)
    chain_metrics = validate_dir / "results" / "chain_metrics.csv"
    system_metrics = validate_dir / "results" / "system_metrics.csv"
    assert_file_exists(chain_metrics)
    assert_file_exists(system_metrics)
    assert_csv_columns_have_values(system_metrics, ("sample_ranking_score", "avg_plddt"))
    assert_csv_columns_have_values(chain_metrics, ("chain_ptm",))
    assert_file_exists(validate_dir / "raw" / "repeat_1" / "normalized" / "artifacts" / "plddt")
    assert_file_exists(validate_dir / "raw" / "repeat_1" / "normalized" / "artifacts" / "pae")
    assert_file_exists(validate_dir / "raw" / "repeat_1" / "normalized" / "artifacts" / "pde")
    return (validate_output,)


@app.cell
def _(mo, validate_output):
    mo.md(f"""
    ### Validate output\n```text\n{validate_output}\n```
    """)
    return


@app.cell
def _(
    assert_csv_columns_have_values,
    assert_file_exists,
    assert_runner_setup_ready,
    build_screen_command,
    env,
    fixtures,
    mo,
    reset_work_dir,
    run_screen,
    runner,
    scoring,
    stream_cli_in_notebook,
    workspace,
):
    mo.stop(
        not run_screen.value,
        mo.md("Enable `Run screen acceptance step` above before starting the expensive screen run."),
    )
    assert_runner_setup_ready(runner, env=env)
    screen_dir = reset_work_dir(workspace / "screen")
    screen_command = build_screen_command(
        runner=runner,
        wrk_dir=screen_dir,
        system_path=fixtures.system_screen_path,
        options_path=fixtures.options_path,
        variable_csv=fixtures.ligand_csv_path,
        scoring_functions=scoring,
    )
    screen_output = stream_cli_in_notebook(screen_command, workspace, env=env)
    summary_csv = screen_dir / "screen_results.csv"
    merged_csv = screen_dir / "screen_results_with_scores.csv"
    assert_file_exists(summary_csv)
    assert_file_exists(merged_csv)
    assert_csv_columns_have_values(
        merged_csv,
        ("system__sample_ranking_score", "system__avg_plddt", "protein_A__chain_ptm"),
    )
    return (screen_output,)


@app.cell
def _(mo, screen_output):
    mo.md(f"""
    ### Screen output\n```text\n{screen_output}\n```
    """)
    return


@app.cell
def _(
    assert_csv_columns_have_values,
    assert_file_exists,
    assert_runner_setup_ready,
    build_oracle_command,
    env,
    fixtures,
    mo,
    oracle_metric_for_runner,
    oracle_scoring_functions_for_runner,
    reset_work_dir,
    run_oracle,
    runner,
    stream_cli_in_notebook,
    workspace,
):
    mo.stop(
        not run_oracle.value,
        mo.md("Enable `Run oracle acceptance step` above before starting the expensive oracle run."),
    )
    assert_runner_setup_ready(runner, env=env)
    oracle_dir = reset_work_dir(workspace / "oracle")
    oracle_metric = oracle_metric_for_runner(runner)
    oracle_command = build_oracle_command(
        runner=runner,
        wrk_dir=oracle_dir,
        system_path=fixtures.system_path,
        options_path=fixtures.options_path,
        input_smiles="CCO",
        output_metric=oracle_metric,
        scoring_functions=oracle_scoring_functions_for_runner(runner),
    )
    oracle_output = stream_cli_in_notebook(oracle_command, workspace, env=env)
    oracle_csv = oracle_dir / "oracle_result.csv"
    assert_file_exists(oracle_csv)
    assert_csv_columns_have_values(oracle_csv, ("value",))
    return (oracle_output,)


@app.cell
def _(mo, oracle_output):
    mo.md(f"""
    ### Oracle output\n```text\n{oracle_output}\n```
    """)
    return


if __name__ == "__main__":
    app.run()
