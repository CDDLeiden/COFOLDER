from __future__ import annotations

import hashlib
from pathlib import Path
import shutil
import tarfile

from cofolder.tools import install_mmseqs


def _archive(temp_dir: Path, *, include_binary: bool = True) -> Path:
    source_root = temp_dir / "archive-source" / "mmseqs"
    source_root.mkdir(parents=True)
    source = source_root / ("mmseqs" if include_binary else "README")
    source.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    archive = temp_dir / "mmseqs.tar.gz"
    with tarfile.open(archive, "w:gz") as bundle:
        bundle.add(source_root.parent, arcname="mmseqs")
    return archive


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_installs_verified_mmseqs_archive(monkeypatch, temp_dir):
    archive = _archive(temp_dir)
    monkeypatch.setattr(
        install_mmseqs,
        "_download",
        lambda _url, destination: shutil.copyfile(archive, destination),
    )
    install_root = temp_dir / "install root"

    result = install_mmseqs.main(
        [
            "--url",
            "https://example.invalid/mmseqs.tar.gz",
            "--sha256",
            _digest(archive),
            "--install-dir",
            str(install_root),
        ]
    )

    installed = install_root / "bin/mmseqs"
    assert result == 0
    assert installed.read_text(encoding="utf-8") == "#!/bin/sh\nexit 0\n"
    assert installed.stat().st_mode & 0o111


def test_rejects_checksum_mismatch_before_install(monkeypatch, temp_dir, capsys):
    archive = _archive(temp_dir)
    monkeypatch.setattr(
        install_mmseqs,
        "_download",
        lambda _url, destination: shutil.copyfile(archive, destination),
    )
    install_root = temp_dir / "install"

    result = install_mmseqs.main(
        [
            "--url",
            "https://example.invalid/mmseqs.tar.gz",
            "--sha256",
            "0" * 64,
            "--install-dir",
            str(install_root),
        ]
    )

    assert result == 1
    assert "sha256 mismatch" in capsys.readouterr().err
    assert not install_root.exists()


def test_rejects_archive_without_mmseqs(monkeypatch, temp_dir, capsys):
    archive = _archive(temp_dir, include_binary=False)
    monkeypatch.setattr(
        install_mmseqs,
        "_download",
        lambda _url, destination: shutil.copyfile(archive, destination),
    )

    result = install_mmseqs.main(
        [
            "--url",
            "https://example.invalid/mmseqs.tar.gz",
            "--sha256",
            _digest(archive),
            "--install-dir",
            str(temp_dir / "install"),
        ]
    )

    assert result == 1
    assert "not found in archive" in capsys.readouterr().err
