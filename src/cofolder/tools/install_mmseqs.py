"""Install a verified vendored MMseqs2 executable for COFOLDER."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
import hashlib
from pathlib import Path, PurePosixPath
import shutil
import stat
import sys
import tarfile
import tempfile
from urllib.request import urlopen


DEFAULT_URL = "https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz"
DEFAULT_SHA256 = "ea05b7c706a166874f56ff0cbdf14791b0de1b7ae3f44e59d2a78383ff8f8d62"


class _Parser(argparse.ArgumentParser):
    """Keep the source shell helper's exit status for invalid arguments."""

    def error(self, message: str) -> None:
        self.print_usage(sys.stderr)
        self.exit(1, f"{self.prog}: error: {message}\n")


def _parser() -> argparse.ArgumentParser:
    parser = _Parser(
        prog="cofolder-tools install-mmseqs",
        description="Install a vendored MMseqs2 binary for COFOLDER."
    )
    parser.add_argument("--url", default=None, help="MMseqs2 tarball URL.")
    parser.add_argument("--sha256", default=None, help="Expected tarball SHA-256.")
    parser.add_argument(
        "--install-dir",
        type=Path,
        default=Path.home() / ".cofolder/vendor/mmseqs",
        help="Installation root (default: ~/.cofolder/vendor/mmseqs).",
    )
    return parser


def _download(url: str, destination: Path) -> None:
    with urlopen(url) as response, destination.open("wb") as output:
        shutil.copyfileobj(response, output)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _extract_mmseqs(archive: Path, destination: Path) -> None:
    with tarfile.open(archive, "r:gz") as bundle:
        candidates = [
            member
            for member in bundle.getmembers()
            if member.isfile() and PurePosixPath(member.name).name == "mmseqs"
        ]
        if not candidates:
            raise RuntimeError("mmseqs executable not found in archive")
        source = bundle.extractfile(candidates[0])
        if source is None:
            raise RuntimeError("mmseqs executable could not be read from archive")
        destination.parent.mkdir(parents=True, exist_ok=True)
        with source, destination.open("wb") as output:
            shutil.copyfileobj(source, output)
    destination.chmod(
        destination.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    if (args.url is None) != (args.sha256 is None):
        parser.error("--url and --sha256 must be provided together")

    url = args.url or DEFAULT_URL
    expected_sha256 = (args.sha256 or DEFAULT_SHA256).lower()
    print(f"[info] source URL: {url}")
    print(f"[info] expected SHA256: {expected_sha256}")

    try:
        with tempfile.TemporaryDirectory(prefix="cofolder-mmseqs-") as temp_dir:
            archive = Path(temp_dir) / "mmseqs.tar.gz"
            print("[info] downloading MMseqs2 archive")
            _download(url, archive)
            print("[info] verifying SHA256")
            actual_sha256 = _sha256(archive)
            if actual_sha256 != expected_sha256:
                raise RuntimeError(
                    "sha256 mismatch\n"
                    f"expected: {expected_sha256}\n"
                    f"actual:   {actual_sha256}"
                )

            print("[info] extracting archive")
            installed_path = args.install_dir.expanduser() / "bin/mmseqs"
            _extract_mmseqs(archive, installed_path)
    except (OSError, RuntimeError, tarfile.TarError) as exc:
        print(f"[error] {exc}", file=sys.stderr)
        return 1

    print(f"[done] installed mmseqs to {installed_path}")
    print(f"[next] export COFOLDER_MMSEQS_BIN='{installed_path}'")
    print(f"[next] or add '{installed_path.parent}' to PATH")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
