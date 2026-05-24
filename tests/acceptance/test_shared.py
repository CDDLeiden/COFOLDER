from __future__ import annotations

import csv
from pathlib import Path

import pytest

from cofolder.acceptance import shared


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


class TestInstallAndRunnerSelection:
    def test_install_command_for_backend(self):
        assert shared.install_command_for_backend("boltz2") == 'pip install "cofolder[acceptance,boltz2]"'

    @pytest.mark.parametrize(
        ("runner", "metric", "groups"),
        [
            ("boltz1", "confidence_score", ["confidence_metrics"]),
            ("boltz2", "affinity_pred_value", ["affinity_metrics"]),
            ("boltz-community", "affinity_pred_value", ["affinity_metrics"]),
        ],
    )
    def test_oracle_runner_settings(self, runner: str, metric: str, groups: list[str]):
        assert shared.oracle_metric_for_runner(runner) == metric
        assert shared.oracle_scoring_functions_for_runner(runner) == groups

    def test_invalid_runner_raises_value_error(self):
        with pytest.raises(ValueError, match="Unsupported backend 'bad-runner'"):
            shared.install_command_for_backend("bad-runner")

    def test_assert_runner_available_uses_runner_check(self, monkeypatch):
        class FakeRunner:
            def check_availability(self):
                return True, None

        monkeypatch.setattr(shared, "get_runner", lambda runner: FakeRunner())

        assert shared.assert_runner_available("boltz2") == "Runner 'boltz2' is available."

    def test_assert_runner_available_surfaces_failure_message(self, monkeypatch):
        class FakeRunner:
            def check_availability(self):
                return False, "install the backend in a clean environment"

        monkeypatch.setattr(shared, "get_runner", lambda runner: FakeRunner())

        with pytest.raises(AssertionError, match="clean boltz1 environment"):
            shared.assert_runner_available("boltz1")


class TestCommandBuilders:
    def test_build_validate_command_uses_cli_contract(self, temp_dir):
        command = shared.build_validate_command(
            runner="boltz2",
            wrk_dir=temp_dir / "validate",
            system_path=temp_dir / "system.yaml",
            options_path=temp_dir / "options.yaml",
            scoring_functions=["confidence_metrics", "affinity_metrics"],
        )

        assert command[:2] == ["cofolder", "validate"]
        assert "-o" in command
        assert "--runner" in command
        assert "--repeats" in command
        assert command[command.index("--runner") + 1] == "boltz2"
        assert command[command.index("--scoring_functions") + 1 :] == [
            "confidence_metrics",
            "affinity_metrics",
        ]

    def test_build_screen_command_contains_wrapper_arguments(self, temp_dir):
        command = shared.build_screen_command(
            runner="boltz1",
            wrk_dir=temp_dir / "screen",
            system_path=temp_dir / "system_screen.yaml",
            options_path=temp_dir / "options.yaml",
            variable_csv=temp_dir / "ligands.csv",
            scoring_functions=["confidence_metrics"],
        )

        assert command[:2] == ["cofolder", "screen"]
        assert "-c" in command
        assert "--col_id" in command
        assert "--variable" in command
        assert "--col_variable" in command
        assert "--merge_data" in command

    def test_build_oracle_command_contains_supported_metric(self, temp_dir):
        command = shared.build_oracle_command(
            runner="boltz-community",
            wrk_dir=temp_dir / "oracle",
            system_path=temp_dir / "system.yaml",
            options_path=temp_dir / "options.yaml",
            input_smiles="CCO",
            output_metric="affinity_pred_value",
            scoring_functions=["affinity_metrics"],
        )

        assert command[:2] == ["cofolder", "oracle"]
        assert command[command.index("--input_smiles") + 1] == "CCO"
        assert command[command.index("--output_metric") + 1] == "affinity_pred_value"
        assert command[command.index("--aggregate") + 1] == "first"


class TestFixtureMaterialization:
    def test_materialize_acceptance_inputs_copies_packaged_files(self, temp_dir):
        inputs = shared.materialize_acceptance_inputs(temp_dir / "fixtures")

        assert inputs.system_path.exists()
        assert inputs.system_screen_path.exists()
        assert inputs.options_path.exists()
        assert inputs.ligand_csv_path.exists()
        assert "diffusion_samples: 1" in inputs.options_path.read_text(encoding="utf-8")

    def test_reset_work_dir_recreates_clean_directory(self, temp_dir):
        work_dir = temp_dir / "workspace"
        work_dir.mkdir()
        stale_file = work_dir / "stale.txt"
        stale_file.write_text("old", encoding="utf-8")

        result = shared.reset_work_dir(work_dir)

        assert result == work_dir
        assert work_dir.exists()
        assert not stale_file.exists()

    def test_run_cli_in_workspace_forwards_cwd(self, monkeypatch, temp_dir):
        captured: dict[str, object] = {}

        class Result:
            returncode = 0
            stdout = "ok"
            stderr = ""

        def fake_run(command, cwd=None, env=None, capture_output=None, text=None, check=None):
            captured["command"] = command
            captured["cwd"] = cwd
            captured["env"] = env
            return Result()

        monkeypatch.setattr(shared.subprocess, "run", fake_run)

        result = shared.run_cli_in_workspace(["cofolder", "validate"], temp_dir)

        assert captured["command"] == ["cofolder", "validate"]
        assert captured["cwd"] == str(temp_dir)
        assert result.stdout == "ok"


class TestCsvAssertions:
    def test_assert_csv_has_columns_and_values(self, temp_dir):
        csv_path = temp_dir / "filled.csv"
        _write_csv(
            csv_path,
            ["name", "affinity_pred_value", "pIC50"],
            [{"name": "row-1", "affinity_pred_value": "1.2", "pIC50": "6.5"}],
        )

        shared.assert_csv_has_columns(csv_path, ["name", "pIC50"])
        shared.assert_csv_columns_have_values(csv_path, ["affinity_pred_value", "pIC50"])

    def test_assert_csv_columns_all_empty(self, temp_dir):
        csv_path = temp_dir / "empty.csv"
        _write_csv(
            csv_path,
            ["name", "affinity_pred_value", "pIC50"],
            [{"name": "row-1", "affinity_pred_value": "", "pIC50": " "}],
        )

        shared.assert_csv_columns_all_empty(csv_path, ["affinity_pred_value", "pIC50"])

    def test_assert_csv_columns_have_values_raises_when_empty(self, temp_dir):
        csv_path = temp_dir / "missing-values.csv"
        _write_csv(
            csv_path,
            ["name", "affinity_pred_value"],
            [{"name": "row-1", "affinity_pred_value": ""}],
        )

        with pytest.raises(AssertionError, match="affinity_pred_value"):
            shared.assert_csv_columns_have_values(csv_path, ["affinity_pred_value"])

    def test_assert_output_contains(self):
        shared.assert_output_contains("hello boltz", "boltz")

        with pytest.raises(AssertionError, match="missing"):
            shared.assert_output_contains("hello boltz", "missing")
