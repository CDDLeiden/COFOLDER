from __future__ import annotations

from cofolder.modules.runners.boltz_runner import BoltzRunner
from cofolder.modules.runners.contracts import RunnerInputCapabilities


class BoltzCommunityRunner(BoltzRunner):
    """Runner for the community-maintained Boltz fork."""

    name = "boltz-community"
    backend_name = "boltz-community"
    backend_distribution = "boltz-community"
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
        available, message = self.check_distribution_available(
            distribution_name="boltz-community",
            missing_message=(
                "The selected 'boltz-community' runner is not installed. Install it with "
                '`pip install "cofolder[boltz-community]"` or '
                '`pip install -e ".[boltz-community]"`.'
            ),
            conflicting_distributions=("boltz",),
            conflict_message=(
                "The selected 'boltz-community' runner cannot be used while a PyPI 'boltz' "
                "installation is present in the same environment because both provide the "
                "`boltz` CLI/module. "
                "Use an environment with only one Boltz-family backend installed."
            ),
        )
        if not available:
            return available, message

        installed_version = self.get_distribution_version("boltz-community")
        if installed_version != "2.10.12":
            return False, (
                "The selected 'boltz-community' runner requires "
                "boltz-community==2.10.12, but boltz-community "
                f"{installed_version} is installed. Install it with "
                '`pip install "cofolder[boltz-community]"` or '
                '`pip install -e ".[boltz-community]"`.'
            )
        return True, None


RUNNER = BoltzCommunityRunner()
