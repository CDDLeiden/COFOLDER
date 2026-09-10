from __future__ import annotations

from cofolder.modules.runners.boltz_runner import BoltzRunner
from cofolder.modules.runners.contracts import RunnerInputCapabilities


class Boltz2Runner(BoltzRunner):
    """Runner for the current Boltz-2 package line."""

    name = "boltz2"
    backend_name = "boltz"
    backend_distribution = "boltz"
    capabilities = {
        "confidence_metrics",
        "affinity_metrics",
        "affinity_metrics_ext",
    }
    model_name = "boltz2"
    input_capabilities = RunnerInputCapabilities(
        entity_types=frozenset({"protein", "ligand", "dna", "rna"}),
        constraint_types=frozenset({"bond", "pocket", "contact"}),
    )

    def check_availability(self) -> tuple[bool, str | None]:
        compatible, message = self._check_python_compatibility(
            runner_name=self.name,
            maximum_exclusive=(3, 13),
        )
        if not compatible:
            return compatible, message

        available, message = self.check_distribution_available(
            distribution_name="boltz",
            missing_message=(
                "The selected 'boltz2' runner is not installed. Install it with "
                '`pip install "cofolder[boltz2]"` or `pip install -e ".[boltz2]"`.'
            ),
            conflicting_distributions=("boltz-community",),
            conflict_message=(
                "The selected 'boltz2' runner cannot be used while 'boltz-community' is installed "
                "in the same environment because both provide the `boltz` CLI/module. "
                "Use an environment with only one Boltz-family backend installed."
            ),
        )
        if not available:
            return available, message

        installed_version = self.get_distribution_version("boltz")
        if installed_version is None:
            return False, (
                "The selected 'boltz2' runner is not installed. Install it with "
                '`pip install "cofolder[boltz2]"` or `pip install -e ".[boltz2]"`.'
            )
        major_version = installed_version.split(".", 1)[0]
        if major_version != "2":
            return False, (
                f"The selected 'boltz2' runner requires the Boltz-2 package line, but boltz "
                f"{installed_version} is installed. Install it with "
                '`pip install "cofolder[boltz2]"` or `pip install -e ".[boltz2]"`.'
            )
        return True, None


RUNNER = Boltz2Runner()
