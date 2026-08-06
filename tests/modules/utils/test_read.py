"""Tests for cofolder.modules.utils.read module."""

from __future__ import annotations

import logging

from cofolder.modules.utils import read


class TestReadJsonLogging:
    def test_read_json_logs_path_at_info_without_dumping_contents(self, temp_dir, caplog):
        json_path = temp_dir / "payload.json"
        json_path.write_text('{"plddt": [1, 2, 3], "score": 0.5}', encoding="utf-8")

        with caplog.at_level(logging.INFO):
            data = read.read_json(json_path)

        assert data == {"plddt": [1, 2, 3], "score": 0.5}
        assert f"{json_path} loaded successfully." in caplog.text
        assert "Contents:" not in caplog.text
        assert "\tplddt:" not in caplog.text
        assert "\tscore:" not in caplog.text

    def test_read_json_logs_contents_at_debug(self, temp_dir, caplog):
        json_path = temp_dir / "payload.json"
        json_path.write_text('{"plddt": [1, 2, 3], "score": 0.5}', encoding="utf-8")

        with caplog.at_level(logging.DEBUG):
            read.read_json(json_path)

        assert f"{json_path} loaded successfully." in caplog.text
        assert f"{json_path} contents:" in caplog.text
        assert "\tplddt: [1, 2, 3]" in caplog.text
        assert "\tscore: 0.5" in caplog.text
