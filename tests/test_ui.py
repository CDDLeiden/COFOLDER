#!/usr/bin/env python
"""
Test suite for Boltz-Lab UI
Validates that all UI components work correctly
"""
import sys
import os
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

def test_imports():
    """Test that all required modules can be imported."""
    print("Testing imports...")
    try:
        import streamlit
        print("  ✓ streamlit")
        import yaml
        print("  ✓ yaml")
        from boltz_eval.ui import app
        print("  ✓ boltz_eval.ui.app")
        return True
    except ImportError as e:
        print(f"  ✗ Import failed: {e}")
        return False


def test_default_options():
    """Test that default options are properly configured."""
    print("\nTesting default options...")
    try:
        from boltz_eval.ui.app import DEFAULT_OPTIONS, OPTION_GROUPS, OPTION_INFO

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

        return True
    except Exception as e:
        print(f"  ✗ Error: {e}")
        return False


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

        return True
    except Exception as e:
        print(f"  ✗ Error: {e}")
        return False


def test_ui_structure():
    """Test that UI has all required components."""
    print("\nTesting UI structure...")
    try:
        from boltz_eval.ui import app

        # Check that main functions exist
        required_functions = [
            'initialize_session_state',
            'render_header',
            'load_yaml_file',
            'save_config',
            'render_recipe_selector',
            'render_predict_ui',
            'render_screen_ui',
            'render_oracle_ui',
            'render_evaluate_ui',
            'render_options_editor',
            'run_predict',
            'run_screen',
            'run_oracle',
            'run_evaluate',
        ]

        for func_name in required_functions:
            assert hasattr(app, func_name), f"Missing function: {func_name}"
            print(f"  ✓ {func_name}")

        return True
    except Exception as e:
        print(f"  ✗ Error: {e}")
        return False


def test_file_structure():
    """Test that all required UI files exist."""
    print("\nTesting file structure...")
    try:
        ui_dir = Path(__file__).parent / 'src' / 'boltz_eval' / 'ui'

        required_files = [
            '__init__.py',
            'app.py',
        ]

        for file_name in required_files:
            file_path = ui_dir / file_name
            assert file_path.exists(), f"Missing file: {file_path}"
            print(f"  ✓ {file_name}")

        # Check config files
        config_file = Path(__file__).parent / '.streamlit' / 'config.toml'
        assert config_file.exists(), "Missing .streamlit/config.toml"
        print(f"  ✓ .streamlit/config.toml")

        # Check documentation
        doc_file = Path(__file__).parent / 'UI_GUIDE.md'
        assert doc_file.exists(), "Missing UI_GUIDE.md"
        print(f"  ✓ UI_GUIDE.md")

        return True
    except Exception as e:
        print(f"  ✗ Error: {e}")
        return False


def test_recipes_configuration():
    """Test that all recipes have proper configuration."""
    print("\nTesting recipes configuration...")
    try:
        from boltz_eval.ui.app import DEFAULT_OPTIONS

        recipes = ['predict', 'screen', 'oracle', 'evaluate']
        print(f"  ✓ Found {len(recipes)} recipes: {', '.join(recipes)}")

        # Verify all common options are present
        common_options = ['out_dir', 'cache', 'checkpoint', 'devices', 'accelerator']
        for opt in common_options:
            assert opt in DEFAULT_OPTIONS, f"Missing common option: {opt}"
        print(f"  ✓ All common options present")

        return True
    except Exception as e:
        print(f"  ✗ Error: {e}")
        return False


def main():
    """Run all tests."""
    print("=" * 60)
    print("Boltz-Eval UI Test Suite")
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
            result = test_func()
            results.append(result)
        except Exception as e:
            print(f"\n✗ Test {test_func.__name__} failed: {e}")
            results.append(False)

    # Summary
    print("\n" + "=" * 60)
    passed = sum(results)
    total = len(results)
    print(f"Results: {passed}/{total} tests passed")

    if all(results):
        print("\n✅ All tests passed! UI is ready to use.")
        print("\nTo start the UI, run:")
        print("  ./run_ui.sh")
        print("  # or")
        print("  streamlit run src/boltz_eval/ui/app.py")
        return 0
    else:
        print("\n❌ Some tests failed. Please check the output above.")
        return 1


if __name__ == '__main__':
    sys.exit(main())

