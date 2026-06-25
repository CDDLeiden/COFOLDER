import marimo

__generated_with = "0.23.8"
app = marimo.App()


@app.cell
def _():
    import sys
    from dataclasses import fields
    from pathlib import Path

    notebook_dir = Path(__file__).resolve().parent
    if str(notebook_dir) not in sys.path:
        sys.path.insert(0, str(notebook_dir))

    import marimo as mo

    from _marimo_helpers import REPO_ROOT, code_block, format_command, read_text
    from cofolder.modules.runners import list_runner_names
    from cofolder.modules.runners.contracts import (
        RunnerExecutionRequest,
        RunnerExecutionResult,
        RunnerMetricOutcome,
        RunnerPreparationResult,
        RunnerRuntime,
    )

    return (
        REPO_ROOT,
        RunnerExecutionRequest,
        RunnerExecutionResult,
        RunnerMetricOutcome,
        RunnerPreparationResult,
        RunnerRuntime,
        code_block,
        fields,
        format_command,
        list_runner_names,
        mo,
        read_text,
    )


@app.cell
def _(REPO_ROOT, list_runner_names):
    runners_dir = REPO_ROOT / "src" / "cofolder" / "modules" / "runners"
    runner_files = sorted(path.name for path in runners_dir.glob("*_runner.py"))
    registered_runners = list_runner_names()
    acceptance_dir = REPO_ROOT / "src" / "cofolder" / "acceptance"
    acceptance_notebooks = sorted(acceptance_dir.glob("*_backend_acceptance.py"))
    return acceptance_dir, acceptance_notebooks, registered_runners, runner_files, runners_dir


@app.cell
def _(
    REPO_ROOT,
    RunnerExecutionRequest,
    RunnerExecutionResult,
    RunnerMetricOutcome,
    RunnerPreparationResult,
    RunnerRuntime,
    acceptance_dir,
    acceptance_notebooks,
    code_block,
    fields,
    format_command,
    mo,
    read_text,
    registered_runners,
    runner_files,
    runners_dir,
):
    runner_files_md = "\n".join(f"- `{name}`" for name in runner_files)
    registered_runners_md = "\n".join(f"- `{name}`" for name in registered_runners)
    contract_lines = "\n".join(
        [
            f"- `RunnerRuntime`: {[field.name for field in fields(RunnerRuntime)]}",
            f"- `RunnerPreparationResult`: {[field.name for field in fields(RunnerPreparationResult)]}",
            f"- `RunnerExecutionRequest`: {[field.name for field in fields(RunnerExecutionRequest)]}",
            f"- `RunnerExecutionResult`: {[field.name for field in fields(RunnerExecutionResult)]}",
            f"- `RunnerMetricOutcome`: {[field.name for field in fields(RunnerMetricOutcome)]}",
        ]
    )
    normalized_tree = "\n".join(
        [
            "raw/repeat_<n>/normalized/",
            "├── system_metrics.csv",
            "├── chain_metrics.csv",
            "├── manifest.json",
            "└── structures/",
            "    └── *.cif | *.mmcif | *.pdb",
        ]
    )
    skeleton = "\n".join(
        [
            "from __future__ import annotations",
            "",
            "from cofolder.modules.runners.base import BaseRunner",
            "from cofolder.modules.runners.contracts import (",
            "    RunnerExecutionRequest,",
            "    RunnerExecutionResult,",
            "    RunnerMetricOutcome,",
            ")",
            "",
            "",
            "class ExampleRunner(BaseRunner):",
            '    name = "example"',
            "",
            "    def run(self, request: RunnerExecutionRequest) -> RunnerExecutionResult:",
            '        normalized_dir = request.repeat_dir / "normalized"',
            "        # runner-owned normalization goes here",
            "        return RunnerExecutionResult(",
            "            runner_name=self.name,",
            "            raw_output_dir=request.repeat_dir,",
            "            normalized_dir=normalized_dir,",
            '            structures_dir=normalized_dir / "structures",',
            '            system_metrics_path=normalized_dir / "system_metrics.csv",',
            '            chain_metrics_path=normalized_dir / "chain_metrics.csv",',
            '            manifest_path=normalized_dir / "manifest.json",',
            "            metric_outcomes={",
            '                "confidence_metrics": RunnerMetricOutcome(state="computed"),',
            "            },",
            "        )",
        ]
    )
    docs_path = REPO_ROOT / "docs" / "tutorials" / "runners.md"
    backend_docs_path = REPO_ROOT / "docs" / "tutorials" / "backend-acceptance.md"
    acceptance_notebooks_md = "\n".join(f"- `{path.name}`" for path in acceptance_notebooks)
    launch_commands = "\n".join(
        format_command(["marimo", "edit", path]) for path in acceptance_notebooks
    )
    mo.md(
        f"""
        # Adding A New Runner

        This notebook is the interactive counterpart to `docs/tutorials/runners.md`.

        Runner modules currently in-tree:

        {runner_files_md}

        Registered runner names:

        {registered_runners_md}

        Contract surface:

        {contract_lines}

        Required normalized bundle layout:

        {code_block(normalized_tree, "text")}

        Reference skeleton:

        {code_block(skeleton, "python")}

        ## Backend Acceptance After Runner Changes

        Use the backend acceptance lane before promoting runner, backend, or shared CLI changes.

        Acceptance notebooks live under:

        `{acceptance_dir}`

        Available backend notebooks:

        {acceptance_notebooks_md}

        Typical launch commands:

        {code_block(launch_commands, "bash")}

        The written contributor guide remains the authoritative source:

        {code_block(read_text(docs_path)[:2200], "markdown")}

        Backend acceptance guide excerpt:

        {code_block(read_text(backend_docs_path)[:1800], "markdown")}
        """
    )
    return


if __name__ == "__main__":
    app.run()
