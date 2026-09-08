from __future__ import annotations

from pathlib import Path

from cofolder.modules.input.command import Command
from cofolder.modules.runners.boltz_runner import BoltzRunner
from cofolder.modules.runners.contracts import RunnerInputCapabilities


class Boltz1Runner(BoltzRunner):
    """Runner for the Boltz-1 package line (boltz==1.0.0)."""

    name = "boltz1"
    model_name = "boltz1"
    capabilities = {
        "confidence_metrics",
    }
    input_capabilities = RunnerInputCapabilities(
        entity_types=frozenset({"protein", "ligand", "dna", "rna"}),
        constraint_types=frozenset({"bond", "pocket"}),
        max_pocket_constraints=1,
        required_pocket_distance=6.0,
    )

    def check_availability(self) -> tuple[bool, str | None]:
        available, message = self.check_distribution_available(
            distribution_name="boltz",
            missing_message=(
                "The selected 'boltz1' runner is not installed. Install it with "
                '`pip install "cofolder[boltz1]"` or `pip install -e ".[boltz1]"`.'
            ),
            conflicting_distributions=("boltz-community",),
            conflict_message=(
                "The selected 'boltz1' runner cannot be used while 'boltz-community' is installed "
                "in the same environment because both provide the `boltz` CLI/module. "
                "Use an environment with only one Boltz-family backend installed."
            ),
        )
        if not available:
            return available, message

        installed_version = self.get_distribution_version("boltz")
        if installed_version != "1.0.0":
            return False, (
                f"The selected 'boltz1' runner requires boltz==1.0.0, but boltz "
                f"{installed_version} is installed. Install it with "
                '`pip install "cofolder[boltz1]"` or `pip install -e ".[boltz1]"`.'
            )
        return True, None

    def load_options(self, options_path: Path) -> Command:
        command = super().load_options(options_path)
        self._set_command_option(command, "model", None)
        return command


RUNNER = Boltz1Runner()
