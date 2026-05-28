from __future__ import annotations

import csv
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Iterable

from cofolder.modules.runners import get_runner


PACKAGE_NAME = "cofolder.acceptance"
DATA_PACKAGE = f"{PACKAGE_NAME}.data"
SUPPORTED_BACKENDS = ("boltz1", "boltz2", "boltz-community", "openfold3")
AFFINITY_COLUMNS = (
    "affinity_pred_value",
    "affinity_probability_binary",
    "pIC50",
    "IC50_M",
    "pIC50_kcal_per_mol",
)
SCREEN_AFFINITY_COLUMNS = tuple(f"ligand_B__{column}" for column in AFFINITY_COLUMNS)
DEFAULT_MANUAL_ROOT = Path.cwd() / ".cofolder-acceptance-runs"
_BACKEND_INSTALL_COMMANDS = {
    "boltz1": 'pip install "cofolder[acceptance,boltz1]"',
    "boltz2": 'pip install "cofolder[acceptance,boltz2]"',
    "boltz-community": 'pip install "cofolder[acceptance,boltz-community]"',
    "openfold3": 'python -m pip install -e ".[acceptance,openfold3]"',
}
_BACKEND_ORACLE_METRICS = {
    "boltz1": "confidence_score",
    "boltz2": "affinity_pred_value",
    "boltz-community": "affinity_pred_value",
    "openfold3": "sample_ranking_score",
}
_BACKEND_ORACLE_SCORING = {
    "boltz1": ["confidence_metrics"],
    "boltz2": ["affinity_metrics"],
    "boltz-community": ["affinity_metrics"],
    "openfold3": ["confidence_metrics"],
}
_BACKEND_ACCEPTANCE_SCORING = {
    "boltz1": ["confidence_metrics", "affinity_metrics", "affinity_metrics_ext"],
    "boltz2": ["confidence_metrics", "affinity_metrics", "affinity_metrics_ext"],
    "boltz-community": ["confidence_metrics", "affinity_metrics", "affinity_metrics_ext"],
    "openfold3": ["confidence_metrics"],
}
_BACKEND_OPTIONS_RESOURCES = {
    "boltz1": "options_acceptance.yaml",
    "boltz2": "options_acceptance.yaml",
    "boltz-community": "options_acceptance.yaml",
    "openfold3": "options_openfold3_acceptance.yaml",
}


@dataclass(frozen=True)
class AcceptanceInputs:
    system_path: Path
    system_screen_path: Path
    options_path: Path
    ligand_csv_path: Path


def install_command_for_backend(runner: str) -> str:
    _validate_runner(runner)
    return _BACKEND_INSTALL_COMMANDS[runner]


def oracle_metric_for_runner(runner: str) -> str:
    _validate_runner(runner)
    return _BACKEND_ORACLE_METRICS[runner]


def oracle_scoring_functions_for_runner(runner: str) -> list[str]:
    _validate_runner(runner)
    return list(_BACKEND_ORACLE_SCORING[runner])


def scoring_functions_for_backend_acceptance(runner: str) -> list[str]:
    _validate_runner(runner)
    return list(_BACKEND_ACCEPTANCE_SCORING[runner])


def options_resource_for_backend(runner: str) -> str:
    _validate_runner(runner)
    return _BACKEND_OPTIONS_RESOURCES[runner]


def materialize_acceptance_inputs(
    target_dir: Path,
    *,
    options_resource: str = "options_acceptance.yaml",
) -> AcceptanceInputs:
    target_dir.mkdir(parents=True, exist_ok=True)
    system_path = _write_package_file("system.yaml", target_dir / "system.yaml")
    system_screen_path = _write_package_file("system_screen.yaml", target_dir / "system_screen.yaml")
    options_path = _write_package_file(options_resource, target_dir / "options_acceptance.yaml")
    ligand_csv_path = _write_package_file("ligand_screen.csv", target_dir / "ligand_screen.csv")
    return AcceptanceInputs(
        system_path=system_path,
        system_screen_path=system_screen_path,
        options_path=options_path,
        ligand_csv_path=ligand_csv_path,
    )


def reset_work_dir(path: Path) -> Path:
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def format_command(command: Iterable[str]) -> str:
    return shlex.join([str(part) for part in command])


def run_cli(
    command: Iterable[str],
    *,
    cwd: Path | None = None,
    check: bool = True,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(
        [str(part) for part in command],
        cwd=str(cwd) if cwd is not None else None,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    if check and completed.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit code {completed.returncode}: {format_command(command)}\n"
            f"STDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
        )
    return completed


def run_cli_in_workspace(
    command: Iterable[str],
    workspace: Path,
    *,
    check: bool = True,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    return run_cli(
        command,
        cwd=workspace,
        check=check,
        env=env,
    )


def assert_runner_available(runner: str) -> str:
    _validate_runner(runner)
    runner_impl = get_runner(runner)
    available, message = runner_impl.check_availability()
    if not available:
        detail = message or f"Runner '{runner}' is not available."
        raise AssertionError(
            f"Expected a clean {runner} environment before running expensive acceptance cells. {detail}"
        )
    return message or f"Runner '{runner}' is available."


def assert_runner_setup_ready(runner: str, *, env: dict[str, str] | None = None) -> str:
    _validate_runner(runner)
    if runner != "openfold3":
        return f"Runner '{runner}' does not require additional setup."

    from cofolder.modules.runners.openfold3_runner import check_openfold3_setup_ready

    ready, message = check_openfold3_setup_ready(env)
    if not ready:
        detail = message or f"Runner '{runner}' setup is incomplete."
        raise AssertionError(
            f"Expected {runner} setup to be ready before running expensive acceptance cells. {detail}"
        )
    return message or f"Runner '{runner}' setup is ready."


def build_validate_command(
    *,
    runner: str,
    wrk_dir: Path,
    system_path: Path,
    options_path: Path,
    scoring_functions: list[str],
) -> list[str]:
    _validate_runner(runner)
    return [
        "cofolder",
        "validate",
        "-s",
        str(system_path),
        "-o",
        str(options_path),
        "-w",
        str(wrk_dir),
        "--runner",
        runner,
        "--repeats",
        "1",
        "--scoring_functions",
        *scoring_functions,
    ]


def build_screen_command(
    *,
    runner: str,
    wrk_dir: Path,
    system_path: Path,
    options_path: Path,
    variable_csv: Path,
    scoring_functions: list[str],
) -> list[str]:
    _validate_runner(runner)
    return [
        "cofolder",
        "screen",
        "-s",
        str(system_path),
        "-o",
        str(options_path),
        "-w",
        str(wrk_dir),
        "--runner",
        runner,
        "--repeats",
        "1",
        "--scoring_functions",
        *scoring_functions,
        "-c",
        str(variable_csv),
        "--col_id",
        "Name",
        "--variable",
        "sequences,1,ligand,smiles",
        "--col_variable",
        "SMILES",
        "--merge_data",
        "pIC50",
    ]


def build_oracle_command(
    *,
    runner: str,
    wrk_dir: Path,
    system_path: Path,
    options_path: Path,
    input_smiles: str,
    output_metric: str,
    scoring_functions: list[str],
) -> list[str]:
    _validate_runner(runner)
    return [
        "cofolder",
        "oracle",
        "-s",
        str(system_path),
        "-o",
        str(options_path),
        "-w",
        str(wrk_dir),
        "--runner",
        runner,
        "--repeats",
        "1",
        "--scoring_functions",
        *scoring_functions,
        "--input_smiles",
        input_smiles,
        "--output_metric",
        output_metric,
        "--aggregate",
        "first",
    ]


def assert_file_exists(path: Path) -> None:
    if not path.exists():
        raise AssertionError(f"Expected file does not exist: {path}")


def assert_csv_has_columns(path: Path, expected_columns: Iterable[str]) -> None:
    rows = _read_csv(path)
    if not rows:
        raise AssertionError(f"Expected CSV with at least one row: {path}")
    missing = [column for column in expected_columns if column not in rows[0]]
    if missing:
        raise AssertionError(f"Missing expected columns in {path}: {missing}")


def assert_csv_columns_all_empty(path: Path, columns: Iterable[str]) -> None:
    rows = _read_csv(path)
    assert_csv_has_columns(path, columns)
    for column in columns:
        if any(_has_value(row.get(column)) for row in rows):
            raise AssertionError(f"Expected column '{column}' to remain empty in {path}")


def assert_csv_columns_have_values(path: Path, columns: Iterable[str]) -> None:
    rows = _read_csv(path)
    assert_csv_has_columns(path, columns)
    for column in columns:
        if not any(_has_value(row.get(column)) for row in rows):
            raise AssertionError(f"Expected column '{column}' to contain at least one value in {path}")


def assert_output_contains(output: str, expected_text: str) -> None:
    if expected_text not in output:
        raise AssertionError(f"Expected output to contain {expected_text!r}")


def combined_output(result: subprocess.CompletedProcess[str]) -> str:
    return "\n".join(part for part in (result.stdout, result.stderr) if part)


def _write_package_file(resource_name: str, destination: Path) -> Path:
    source = resources.files(DATA_PACKAGE).joinpath(resource_name)
    destination.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")
    return destination


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _has_value(value: str | None) -> bool:
    if value is None:
        return False
    return str(value).strip() != ""


def _validate_runner(runner: str) -> None:
    if runner not in SUPPORTED_BACKENDS:
        supported = ", ".join(SUPPORTED_BACKENDS)
        raise ValueError(f"Unsupported backend '{runner}'. Expected one of: {supported}")
