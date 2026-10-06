"""Credential-redaction tests for backend command reporting."""

import logging
import subprocess

import pytest

from cofolder.modules.contracts import (
    PUBLIC_SCHEMA_VERSION,
    FailureStage,
    OutputIdentity,
    PublicManifest,
    PublicOutputBundle,
    WorkflowKind,
    failure_from_exception,
    write_public_bundle,
)
from cofolder.modules.runners._command_reporting import REDACTED, redact_sensitive_data
from cofolder.modules.runners.boltz_runner import run_boltz
from cofolder.modules.runners.openfold3_runner import OpenFold3Runner, run_openfold3


class _Output:
    def __init__(self, lines):
        self._lines = lines

    def __iter__(self):
        return iter(self._lines)

    def close(self):
        return None


class _Process:
    def __init__(self, lines, returncode):
        self.stdout = _Output(lines)
        self._returncode = returncode

    def wait(self):
        return self._returncode


def test_recursive_diagnostic_redaction_preserves_non_secret_header():
    assert redact_sensitive_data(
        {
            "runner": {
                "msa_server_password": "password-value",
                "api_key_header": "X-API-Key",
                "nested": [{"access-token": "token-value"}],
            }
        }
    ) == {
        "runner": {
            "msa_server_password": REDACTED,
            "api_key_header": "X-API-Key",
            "nested": [{"access-token": REDACTED}],
        }
    }


@pytest.mark.parametrize("returncode", [0, 9])
def test_boltz_uses_real_secret_but_reports_only_redacted_diagnostics(
    monkeypatch, caplog, tmp_path, returncode
):
    secret = "boltz-dummy-credential"
    actual = {}

    def _popen(argv, **kwargs):
        actual["argv"] = argv
        return _Process([f"backend echoed {secret}\n"], returncode)

    monkeypatch.setattr(
        "cofolder.modules.runners.boltz_runner.subprocess.Popen", _popen
    )
    caplog.set_level(logging.INFO)
    argv = ["boltz", "predict", "system.yaml", "--api_key_value", secret]

    if returncode:
        with pytest.raises(subprocess.CalledProcessError) as caught:
            run_boltz(argv)
        result = caught.value
        failure = failure_from_exception(
            result,
            identity=OutputIdentity(
                workflow=WorkflowKind.VALIDATE,
                run_id="run-1",
                system_id="system-1",
                runner_id="boltz2",
            ),
            stage=FailureStage.BACKEND_EXECUTION,
            error_code="runner_backend_execution_failed",
        )
        output = write_public_bundle(
            PublicOutputBundle(
                manifest=PublicManifest(
                    schema_version=PUBLIC_SCHEMA_VERSION,
                    identity=failure.envelope.identity,
                    status="failed",
                ),
                records=(failure,),
            ),
            tmp_path / "public",
        )
        persisted = output.records_path.read_text() + output.failures_path.read_text()
        assert secret not in persisted
        assert result.returncode == 9
    else:
        result = run_boltz(argv)
        assert result.returncode == 0

    assert actual["argv"] == argv
    assert secret not in str(result.args)
    assert secret not in str(result.stdout)
    assert secret not in caplog.text
    assert REDACTED in str(result.args)
    assert REDACTED in str(result.stdout)


@pytest.mark.parametrize("fails", [False, True])
def test_openfold3_uses_real_secret_but_sanitizes_results_and_errors(
    monkeypatch, caplog, tmp_path, fails
):
    secret = "openfold-dummy-token"
    options_path = tmp_path / "options.yaml"
    options_path.write_text(
        "version: 1\n"
        "runtime:\n"
        "  extra_args:\n"
        f"    - --access-token={secret}\n"
        "runner: {}\n",
        encoding="utf-8",
    )
    options = OpenFold3Runner().load_options(options_path)
    actual = {}

    def _run(argv, **kwargs):
        actual["argv"] = argv
        if fails:
            raise subprocess.CalledProcessError(
                4, argv, output=f"stdout {secret}", stderr=f"stderr {secret}"
            )
        return subprocess.CompletedProcess(
            argv, 0, stdout=f"stdout {secret}", stderr=f"stderr {secret}"
        )

    monkeypatch.setattr(
        "cofolder.modules.runners.openfold3_runner.subprocess.run", _run
    )
    caplog.set_level(logging.INFO)

    def call():
        return run_openfold3(
            query_json_path=tmp_path / "query.json",
            runner_yaml_path=tmp_path / "runner.yaml",
            output_dir=tmp_path / "out",
            diffusion_samples=1,
            options=options,
        )
    if fails:
        with pytest.raises(subprocess.CalledProcessError) as caught:
            call()
        result = caught.value
        assert result.returncode == 4
        assert secret not in str(result.stderr)
    else:
        result = call()
        assert result.returncode == 0

    assert any(secret in token for token in actual["argv"])
    assert secret not in str(result.cmd if fails else result.args)
    assert secret not in str(result.output if fails else result.stdout)
    assert secret not in caplog.text
    assert REDACTED in str(result.cmd if fails else result.args)
