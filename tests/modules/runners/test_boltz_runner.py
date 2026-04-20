"""Tests for Boltz runner timing parsing."""

from cofolder.modules.runners.boltz_runner import _parse_boltz_stage_timings


def test_parse_boltz_stage_timings_with_msa_and_affinity():
    timed_lines = [
        (0.5, "Calling MSA server for target system with 1 sequences"),
        (1.0, "SUBMIT:   0%|          | 0/150 [elapsed: 00:00 remaining: ?]"),
        (4.0, "COMPLETE: 100%|██████████| 150/150 [elapsed: 00:03 remaining: 00:00]"),
        (7.0, "Running affinity prediction for 1 input."),
    ]

    out = _parse_boltz_stage_timings(timed_lines, total_elapsed=15.0)

    assert out["boltz.msa"] == 3.5
    assert out["boltz.affinity_prediction"] == 8.0


def test_parse_boltz_stage_timings_without_msa_or_affinity():
    timed_lines = [
        (0.2, "Checking input data."),
        (3.0, "Predicting structure for 1 input."),
    ]

    out = _parse_boltz_stage_timings(timed_lines, total_elapsed=10.0)

    assert out == {}


def test_parse_boltz_stage_timings_skips_affinity_when_outputs_already_exist():
    timed_lines = [
        (0.4, "Predicting property: affinity"),
        (0.5, "Found some existing affinity predictions (1), skipping and running only the missing ones, if any. If you wish to override these existing affinity predictions, please set the --override flag."),
        (0.6, "Found existing affinity predictions for all inputs, skipping."),
    ]

    out = _parse_boltz_stage_timings(timed_lines, total_elapsed=5.0)

    assert out == {}


def test_parse_boltz_stage_timings_ignores_incomplete_msa_window():
    timed_lines = [
        (0.5, "Calling MSA server for target system with 1 sequences"),
        (2.0, "SUBMIT:   0%|          | 0/150 [elapsed: 00:00 remaining: ?]"),
        (4.0, "Running affinity prediction for 1 input."),
    ]

    out = _parse_boltz_stage_timings(timed_lines, total_elapsed=9.0)

    assert "boltz.msa" not in out
    assert out["boltz.affinity_prediction"] == 5.0
