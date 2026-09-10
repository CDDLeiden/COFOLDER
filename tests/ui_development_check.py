#!/usr/bin/env python
"""Optional checks for the unsupported, obsolete source-tree UI prototype.

This file intentionally does not match pytest's default test filename pattern and
is not part of the supported release test suite. Its checks require Streamlit and
do not establish compatibility with current COFOLDER recipe contracts.
"""
import sys
import os
import importlib
from pathlib import Path
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent

# Add repo src to path
sys.path.insert(0, str(REPO_ROOT / 'src'))

def test_imports():
    """Test that all required modules can be imported."""
    print("Testing imports...")
    try:
        for module_name in ("streamlit", "yaml", "cofolder.ui.app"):
            importlib.import_module(module_name)
            print(f"  ✓ {module_name}")
    except ImportError as e:
        pytest.fail(f"Import failed: {e}")


def test_default_options():
    """Test that default options are properly configured."""
    print("\nTesting default options...")
    try:
        from cofolder.ui.app import DEFAULT_OPTIONS, OPTION_GROUPS, OPTION_INFO

        # Check that DEFAULT_OPTIONS is a dict
        assert isinstance(DEFAULT_OPTIONS, dict), "DEFAULT_OPTIONS should be a dict"
        print(f"  ✓ DEFAULT_OPTIONS has {len(DEFAULT_OPTIONS)} options")

        # Check option groups
        assert isinstance(OPTION_GROUPS, dict), "OPTION_GROUPS should be a dict"
        total_options = sum(len(opts) for opts in OPTION_GROUPS.values())
        print(f"  ✓ OPTION_GROUPS defines {len(OPTION_GROUPS)} groups with {total_options} options")

        # Check option info
        assert isinstance(OPTION_INFO, dict), "OPTION_INFO should be a dict"
        print(f"  ✓ OPTION_INFO has descriptions for {len(OPTION_INFO)} options")

        # Verify all options in groups exist in defaults
        for group, options in OPTION_GROUPS.items():
            for opt in options:
                if opt not in DEFAULT_OPTIONS:
                    print(f"  ⚠ {opt} in {group} but not in DEFAULT_OPTIONS")

    except Exception as e:
        pytest.fail(f"Error: {e}")


def test_yaml_operations():
    """Test YAML loading and saving functionality."""
    print("\nTesting YAML operations...")
    try:
        import yaml
        import tempfile

        # Create test config
        test_config = {
            'wrapper': [
                {'run_dir': './test'},
                {'system': './test.yaml'},
            ],
            'options': [{
                'out_dir': './output',
                'cache': '~/.boltz',
                'devices': 1,
            }]
        }

        # Test saving
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(test_config, f)
            temp_file = f.name

        # Test loading
        with open(temp_file, 'r') as f:
            loaded_config = yaml.safe_load(f)

        assert loaded_config == test_config, "Loaded config doesn't match saved config"
        print("  ✓ YAML save/load works correctly")

        # Cleanup
        os.unlink(temp_file)

    except Exception as e:
        pytest.fail(f"Error: {e}")


def test_ui_structure():
    """Test that UI has all required components."""
    print("\nTesting UI structure...")
    try:
        from cofolder.ui import app

        # Check that main functions exist
        required_functions = [
            'initialize_session_state',
            'render_header',
            'load_yaml_file',
            'save_config',
            'render_recipe_selector',
            'render_validate_ui',
            'render_screen_ui',
            'render_oracle_ui',
            'render_options_editor',
            'run_validate',
            'run_screen',
            'run_oracle',
        ]

        for func_name in required_functions:
            assert hasattr(app, func_name), f"Missing function: {func_name}"
            print(f"  ✓ {func_name}")

    except Exception as e:
        pytest.fail(f"Error: {e}")


def test_file_structure():
    """Test that all required UI files exist."""
    print("\nTesting file structure...")
    try:
        ui_dir = REPO_ROOT / 'src' / 'cofolder' / 'ui'

        required_files = [
            '__init__.py',
            'app.py',
        ]

        for file_name in required_files:
            file_path = ui_dir / file_name
            assert file_path.exists(), f"Missing file: {file_path}"
            print(f"  ✓ {file_name}")

        # Check config files
        config_file = REPO_ROOT / '.streamlit' / 'config.toml'
        assert config_file.exists(), "Missing .streamlit/config.toml"
        print("  ✓ .streamlit/config.toml")

        # Check launcher script
        launcher_file = REPO_ROOT / 'run_ui.sh'
        assert launcher_file.exists(), "Missing run_ui.sh"
        print("  ✓ run_ui.sh")

    except Exception as e:
        pytest.fail(f"Error: {e}")


def test_recipes_configuration():
    """Test that all recipes have proper configuration."""
    print("\nTesting recipes configuration...")
    try:
        from cofolder.ui.app import DEFAULT_OPTIONS

        recipes = ['validate', 'screen', 'oracle']
        print(f"  ✓ Found {len(recipes)} recipes: {', '.join(recipes)}")

        # Verify all common options are present
        common_options = ['out_dir', 'cache', 'checkpoint', 'devices', 'accelerator']
        for opt in common_options:
            assert opt in DEFAULT_OPTIONS, f"Missing common option: {opt}"
        print("  ✓ All common options present")

    except Exception as e:
        pytest.fail(f"Error: {e}")


def main():
    """Run the optional development checks."""
    print("=" * 60)
    print("UNSUPPORTED COFOLDER UI DEVELOPMENT CHECKS")
    print("These checks do not establish compatibility with current recipe contracts.")
    print("=" * 60)

    tests = [
        test_imports,
        test_default_options,
        test_yaml_operations,
        test_ui_structure,
        test_file_structure,
        test_recipes_configuration,
    ]

    results = []
    for test_func in tests:
        try:
            test_func()
            results.append(True)
        except BaseException as e:
            print(f"\n✗ Test {test_func.__name__} failed: {e}")
            results.append(False)

    # Summary
    print("\n" + "=" * 60)
    passed = sum(results)
    total = len(results)
    print(f"Results: {passed}/{total} tests passed")

    if all(results):
        print("\nAll optional development checks passed.")
        print("The UI remains unsupported and outside the COFOLDER release surface.")
        return 0
    else:
        print("\n❌ Some tests failed. Please check the output above.")
        return 1


if __name__ == '__main__':
    sys.exit(main())
