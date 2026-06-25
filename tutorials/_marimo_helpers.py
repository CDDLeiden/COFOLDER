from __future__ import annotations

from pathlib import Path
import shlex
import subprocess
import sys
import tempfile


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
EXAMPLES_DIR = REPO_ROOT / "examples"
DOCS_DIR = REPO_ROOT / "docs"
TUTORIALS_DIR = REPO_ROOT / "tutorials"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))


def relative_to_repo(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def code_block(text: str, language: str = "") -> str:
    stripped = text.rstrip()
    return f"```{language}\n{stripped}\n```"


def format_command(command: list[str | Path]) -> str:
    return shlex.join(str(part) for part in command)


def run_command(
    command: list[str | Path],
    *,
    cwd: Path | None = None,
) -> str:
    resolved = [str(part) for part in command]
    result = subprocess.run(
        resolved,
        cwd=None if cwd is None else str(cwd),
        capture_output=True,
        text=True,
        check=False,
    )
    output = "\n".join(
        piece for piece in [result.stdout.strip(), result.stderr.strip()] if piece
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit code {result.returncode}: {format_command(resolved)}\n"
            f"{output}"
        )
    return output


def make_workspace(prefix: str) -> Path:
    return Path(tempfile.mkdtemp(prefix=prefix))
