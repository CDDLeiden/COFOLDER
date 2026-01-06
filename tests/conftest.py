"""Pytest configuration and fixtures for boltz-lab tests."""
import tempfile
from pathlib import Path

import pytest
import yaml


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_yaml_file(temp_dir):
    """Create a sample YAML file for testing."""
    yaml_path = temp_dir / "test.yaml"
    data = {"key1": "value1", "key2": {"nested": "value2"}}
    with open(yaml_path, "w") as f:
        yaml.dump(data, f)
    return yaml_path


@pytest.fixture
def sample_system_yaml(temp_dir):
    """Create a sample system YAML file."""
    system_data = {
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "fasta": "MKRAAT"
                }
            },
            {
                "ligand": {
                    "smiles": "CCO",
                    "ccd": "ETH"
                }
            }
        ]
    }
    yaml_path = temp_dir / "system.yaml"
    with open(yaml_path, "w") as f:
        yaml.dump(system_data, f)
    return yaml_path


@pytest.fixture
def sample_options_yaml(temp_dir):
    """Create a sample Boltz options YAML file."""
    options_data = {
        "options": [
            {"cache": "~/.boltz"},
            {"recycling_steps": 3},
            {"diffusion_samples": 1}
        ]
    }
    yaml_path = temp_dir / "options.yaml"
    with open(yaml_path, "w") as f:
        yaml.dump(options_data, f)
    return yaml_path


@pytest.fixture
def sample_csv_file(temp_dir):
    """Create a sample CSV file for testing."""
    csv_path = temp_dir / "test.csv"
    with open(csv_path, "w") as f:
        f.write("compound_id,smiles,mw\n")
        f.write("CMPD001,CCO,46.07\n")
        f.write("CMPD002,CC(C)O,60.10\n")
    return csv_path
