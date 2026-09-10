"""Secret-safe diagnostic representations of external commands."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

REDACTED = "[REDACTED]"


def _normalized_key(value: object) -> str:
    return str(value).lstrip("-").lower().replace("-", "_")


def is_sensitive_key(value: object) -> bool:
    """Return whether an option name conventionally carries a credential."""

    key = _normalized_key(value)
    if key == "api_key_header":
        return False
    return (
        key == "api_key_value"
        or key.endswith("password")
        or key.endswith("secret")
        or key.endswith("token")
    )


def redact_sensitive_data(value: Any) -> Any:
    """Recursively redact credential-bearing values in diagnostic data."""

    if isinstance(value, Mapping):
        return {
            key: REDACTED if is_sensitive_key(key) else redact_sensitive_data(child)
            for key, child in value.items()
        }
    if isinstance(value, tuple):
        return tuple(redact_sensitive_data(child) for child in value)
    if isinstance(value, list):
        return [redact_sensitive_data(child) for child in value]
    return value


@dataclass(frozen=True, slots=True)
class CommandReport:
    """A redacted argv plus the secret values needed to sanitize output."""

    argv: tuple[str, ...]
    secret_values: tuple[str, ...]

    def redact_text(self, value: object | None) -> str | None:
        if value is None:
            return None
        result = str(value)
        for secret in sorted(self.secret_values, key=len, reverse=True):
            if secret:
                result = result.replace(secret, REDACTED)
        return result


def command_report(argv: Sequence[object]) -> CommandReport:
    """Build a display-only argv while retaining no credential values in it."""

    tokens = [str(token) for token in argv]
    display: list[str] = []
    secrets: list[str] = []
    index = 0
    while index < len(tokens):
        token = tokens[index]
        if token.startswith("--") and "=" in token:
            name, value = token.split("=", 1)
            if is_sensitive_key(name):
                display.append(f"{name}={REDACTED}")
                if value:
                    secrets.append(value)
            else:
                display.append(token)
            index += 1
            continue

        display.append(token)
        if token.startswith("--") and is_sensitive_key(token) and index + 1 < len(tokens):
            secret = tokens[index + 1]
            display.append(REDACTED)
            if secret:
                secrets.append(secret)
            index += 2
            continue
        index += 1

    return CommandReport(tuple(display), tuple(dict.fromkeys(secrets)))
