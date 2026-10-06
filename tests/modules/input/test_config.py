from __future__ import annotations

from pathlib import Path

import pytest

from cofolder.modules.input import OptionsValidationError, YamlLoadError
from cofolder.modules.input.config import (
    BOLTZ_OPTIONS_SCHEMA,
    OPENFOLD3_OPTIONS_SCHEMA,
    load_runner_options,
    load_yaml_document,
)


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_yaml_syntax_error_includes_source_location(temp_dir):
    path = _write(temp_dir / "bad.yaml", "version: [1\n")

    with pytest.raises(YamlLoadError) as caught:
        load_yaml_document(path)

    assert caught.value.source_path == path
    assert caught.value.line is not None and caught.value.line >= 1
    assert caught.value.column is not None
    assert caught.value.error_code == "yaml_load_failed"


@pytest.mark.parametrize(
    ("text", "field_path"),
    [
        ("runner: {}\n", ("version",)),
        ("version: 1\n", ("runner",)),
        ("version: 1\nrunner: {}\nunknown: true\n", ("unknown",)),
        ("version: 1\nruntime:\n  diffusion_samples: many\nrunner: {}\n", ("runtime", "diffusion_samples")),
        ("version: 1\nrunner:\n  accelerator: quantum\n", ("runner", "accelerator")),
    ],
)
def test_boltz_options_reject_invalid_schema(temp_dir, text, field_path):
    path = _write(temp_dir / "options.yaml", text)

    with pytest.raises(OptionsValidationError) as caught:
        load_runner_options(path, schema=BOLTZ_OPTIONS_SCHEMA)

    assert caught.value.source_path == path
    assert caught.value.field_path == field_path
    assert caught.value.error_code == "options_validation_failed"


def test_boltz_options_load_namespaced_values(temp_dir):
    path = _write(
        temp_dir / "options.yaml",
        "version: 1\nruntime:\n  cache_path: ./cache\n  diffusion_samples: 2\n"
        "runner:\n  accelerator: cpu\n  use_potentials: true\n",
    )

    options = load_runner_options(path, schema=BOLTZ_OPTIONS_SCHEMA)

    assert options.version == 1
    assert options.runtime.cache_path == Path("cache")
    assert options.runtime.diffusion_samples == 2
    assert options.runner == {"accelerator": "cpu", "use_potentials": True}


def test_openfold_options_reject_unknown_runner_key(temp_dir):
    path = _write(
        temp_dir / "options.yaml",
        "version: 1\nruntime: {}\nrunner:\n  arbitrary_passthrough: true\n",
    )

    with pytest.raises(OptionsValidationError, match="arbitrary_passthrough"):
        load_runner_options(path, schema=OPENFOLD3_OPTIONS_SCHEMA)


@pytest.mark.parametrize(
    ("fragment", "field_path"),
    [
        (
            "experiment_settings:\n    seeds: [42]",
            ("runner", "experiment_settings", "seeds"),
        ),
        (
            "experiment_settings:\n    use_msa_server: true",
            ("runner", "experiment_settings", "use_msa_server"),
        ),
        (
            "output_writer_settings:\n    unknown_output: true",
            ("runner", "output_writer_settings", "unknown_output"),
        ),
    ],
)
def test_openfold_options_reject_unknown_or_workflow_owned_nested_keys(
    temp_dir, fragment, field_path
):
    path = _write(
        temp_dir / "options.yaml",
        f"version: 1\nruntime: {{}}\nrunner:\n  {fragment}\n",
    )

    with pytest.raises(OptionsValidationError) as caught:
        load_runner_options(path, schema=OPENFOLD3_OPTIONS_SCHEMA)

    assert caught.value.field_path == field_path


def test_openfold_options_accept_explicit_native_configuration_tree(temp_dir):
    path = _write(
        temp_dir / "options.yaml",
        """version: 1
runtime:
  diffusion_samples: 2
runner:
  experiment_settings:
    use_templates: false
  pl_trainer_args:
    accelerator: gpu
    devices: 1
  model_update:
    presets: [predict, low_mem]
    custom:
      settings:
        memory:
          eval:
            use_cueq_triangle_kernels: false
  output_writer_settings:
    structure_format: cif
    full_confidence_output_dtype: float16
""",
    )

    options = load_runner_options(path, schema=OPENFOLD3_OPTIONS_SCHEMA)

    assert options.runner["pl_trainer_args"]["devices"] == 1
    assert options.runner["model_update"]["presets"] == ["predict", "low_mem"]
