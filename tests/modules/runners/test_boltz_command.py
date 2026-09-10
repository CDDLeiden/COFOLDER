"""Regression tests for typed Boltz argument construction."""

import pytest
import yaml

from cofolder.modules.input.config import RunnerOptions
from cofolder.modules.runners._boltz_command import build_boltz_command
from cofolder.modules.runners.boltz1_runner import Boltz1Runner
from cofolder.modules.runners.boltz2_runner import Boltz2Runner
from cofolder.modules.runners.boltz_community_runner import BoltzCommunityRunner


@pytest.mark.parametrize(
    ("runner", "expected_model"),
    [
        (Boltz1Runner(), []),
        (Boltz2Runner(), ["--model", "boltz2"]),
        (BoltzCommunityRunner(), ["--model", "boltz2"]),
    ],
)
def test_load_options_builds_complete_typed_argv(
    runner, expected_model, temp_dir
):
    options_path = temp_dir / f"{runner.name}.yaml"
    # Deliberately use non-schema order to prove source key order is irrelevant.
    options_path.write_text(
        yaml.safe_dump(
            {
                "runner": {
                    "override": False,
                    "write_full_pae": True,
                    "step_scale": 0.0,
                    "recycling_steps": 0,
                    "devices": 1,
                },
                "runtime": {
                    "diffusion_samples": 2,
                    "cache_path": "/typed/cache",
                },
                "version": 1,
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )

    options = runner.load_options(options_path)

    assert isinstance(options, RunnerOptions)
    assert build_boltz_command(
        options=options,
        system_path="system.yaml",
        output_dir="output",
        seed=7,
        model_name=runner.command_model_name,
        use_msa_server=True,
    ) == [
        "boltz",
        "predict",
        "system.yaml",
        "--out_dir",
        "output",
        "--seed",
        "7",
        "--use_msa_server",
        "--cache",
        "/typed/cache",
        "--diffusion_samples",
        "2",
        "--devices",
        "1",
        "--recycling_steps",
        "0",
        "--step_scale",
        "0.0",
        "--write_full_pae",
        *expected_model,
    ]


def test_typed_argv_uses_default_sample_count_and_omits_absent_fields(temp_dir):
    options_path = temp_dir / "options.yaml"
    options_path.write_text("version: 1\nruntime: {}\nrunner: {}\n", encoding="utf-8")
    runner = Boltz2Runner()

    argv = build_boltz_command(
        options=runner.load_options(options_path),
        system_path="system.yaml",
        output_dir="output",
        seed=0,
        model_name=runner.model_name,
        use_msa_server=False,
    )

    assert argv == [
        "boltz",
        "predict",
        "system.yaml",
        "--out_dir",
        "output",
        "--seed",
        "0",
        "--diffusion_samples",
        "1",
        "--model",
        "boltz2",
    ]
