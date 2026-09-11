import marimo

__generated_with = "0.23.8"
app = marimo.App()


@app.cell
def _():
    import subprocess
    import tempfile
    from pathlib import Path

    import marimo as mo

    from cofolder.acceptance.shared import (
        AFFINITY_COLUMNS,
        SCREEN_MANUSCRIPT_COLUMNS,
        assert_csv_columns_have_values,
        assert_file_exists,
        assert_runner_available,
        assert_screen_msa_reused_once,
        build_oracle_command,
        build_screen_command,
        build_validate_command,
        format_command,
        install_command_for_backend,
        materialize_acceptance_inputs,
        oracle_metric_for_runner,
        oracle_scoring_functions_for_runner,
        reset_work_dir,
        scoring_functions_for_backend_acceptance,
    )

    def stream_cli_in_notebook(command: list[str], workspace: Path) -> str:
        print(f"$ {format_command(command)}", flush=True)
        process = subprocess.Popen(
            command,
            cwd=str(workspace),
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
        AFFINITY_COLUMNS,
        Path,
        SCREEN_MANUSCRIPT_COLUMNS,
        assert_csv_columns_have_values,
        assert_file_exists,
        assert_runner_available,
        assert_screen_msa_reused_once,
        build_oracle_command,
        build_screen_command,
        build_validate_command,
        format_command,
        install_command_for_backend,
        materialize_acceptance_inputs,
        mo,
        oracle_metric_for_runner,
        oracle_scoring_functions_for_runner,
        reset_work_dir,
        scoring_functions_for_backend_acceptance,
        stream_cli_in_notebook,
        tempfile,
    )


@app.cell
def _(mo):
    mo.md("""
    # Boltz2 backend acceptance

    Install command for a fresh environment:

    ```bash
    pip install "cofolder[acceptance,analysis,boltz2]"
    ```

    This manual lane is intentionally expensive and is **not** part of routine `pytest` or CI-default checks.
    Run each expensive command deliberately by enabling its checkbox: `validate`, then `screen`, then `oracle`.

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
    tempfile,
):
    runner = "boltz2"
    workspace = Path(tempfile.mkdtemp(prefix="cofolder-boltz2-acceptance-"))
    fixtures = materialize_acceptance_inputs(workspace / "fixtures")
    cache_dir = workspace / ".boltz" / "cache"
    cache_dir.parent.mkdir(parents=True, exist_ok=True)
    options_text = fixtures.options_path.read_text(encoding="utf-8")
    cache_placeholder = "./cache/.boltz"
    if cache_placeholder not in options_text:
        raise AssertionError(
            f"Expected cache placeholder {cache_placeholder!r} in {fixtures.options_path}, "
            "but it was not found. Update the acceptance notebook cache rewrite."
        )
    options_text = options_text.replace(cache_placeholder, str(cache_dir), 1)
    fixtures.options_path.write_text(options_text, encoding="utf-8")
    install_command = install_command_for_backend(runner)
    availability_message = assert_runner_available(runner)
    return (
        availability_message,
        cache_dir,
        fixtures,
        install_command,
        runner,
        workspace,
    )


@app.cell
def _(availability_message, cache_dir, install_command, mo):
    mo.md(f"""
    **Verified install command:** `{install_command}`

    **Availability check:** `{availability_message}`

    **Per-run cache path:** `{cache_dir}`

    The packaged options file was rewritten for this run so the backend cache stays inside the temporary workspace.
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

        **Options cache path:** `{workspace / ".boltz" / "cache"}`

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
    AFFINITY_COLUMNS,
    assert_csv_columns_have_values,
    assert_file_exists,
    build_validate_command,
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
    validate_dir = reset_work_dir(workspace / "validate")
    validate_command = build_validate_command(
        runner=runner,
        wrk_dir=validate_dir,
        system_path=fixtures.system_path,
        options_path=fixtures.options_path,
        scoring_functions=scoring,
    )
    validate_output = stream_cli_in_notebook(validate_command, workspace)
    chain_metrics = validate_dir / "results" / "metrics.csv"
    system_metrics = chain_metrics
    assert_file_exists(chain_metrics)
    assert_file_exists(system_metrics)
    assert_csv_columns_have_values(chain_metrics, AFFINITY_COLUMNS)
    return (validate_output,)


@app.cell
def _(mo, validate_output):
    mo.md(f"""
    ### Validate output\n```text\n{validate_output}\n```
    """)
    return


@app.cell
def _(
    SCREEN_MANUSCRIPT_COLUMNS,
    assert_csv_columns_have_values,
    assert_file_exists,
    assert_screen_msa_reused_once,
    build_screen_command,
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
    screen_dir = reset_work_dir(workspace / "screen")
    screen_command = build_screen_command(
        runner=runner,
        wrk_dir=screen_dir,
        system_path=fixtures.system_screen_path,
        options_path=fixtures.options_path,
        library=fixtures.ligand_csv_path,
        scoring_functions=scoring,
        protein_training_data_path=fixtures.protein_training_data_path,
        ligand_training_data_path=fixtures.ligand_training_data_path,
        pocket_coverage_reference="F1",
    )
    screen_command.append("--cluster_ifps")
    screen_output = stream_cli_in_notebook(screen_command, workspace)
    summary_csv = screen_dir / "results" / "metrics.csv"
    merged_csv = summary_csv
    cluster_summary_csv = screen_dir / "results" / "ifp_cluster_summary.csv"
    assert_file_exists(summary_csv)
    assert_file_exists(merged_csv)
    assert_file_exists(cluster_summary_csv)
    assert_csv_columns_have_values(merged_csv, SCREEN_MANUSCRIPT_COLUMNS)
    assert_csv_columns_have_values(
        summary_csv,
        ("ifp_cluster_id", "ifp_cluster_status"),
    )
    msa_reuse_message = assert_screen_msa_reused_once(screen_dir, runner)
    return msa_reuse_message, screen_output


@app.cell
def _(mo, msa_reuse_message, screen_output):
    mo.md(f"""
    **MSA reuse check:** `{msa_reuse_message}`

    ### Screen output\n```text\n{screen_output}\n```
    """)
    return


@app.cell
def _(
    assert_csv_columns_have_values,
    assert_file_exists,
    build_oracle_command,
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
    oracle_output = stream_cli_in_notebook(oracle_command, workspace)
    oracle_csv = oracle_dir / "results" / "metrics.csv"
    assert_file_exists(oracle_csv)
    assert_csv_columns_have_values(oracle_csv, ("oracle_score",))
    return (oracle_output,)


@app.cell
def _(mo, oracle_output):
    mo.md(f"""
    ### Oracle output\n```text\n{oracle_output}\n```
    """)
    return


if __name__ == "__main__":
    app.run()
