"""Tests for cofolder.recipes.validate module."""

from cofolder.recipes.validate import Validate


class TestValidateInit:
    """Tests for Validate initialization."""

    def test_init_basic(self, sample_system_yaml, sample_options_yaml, temp_dir):
        """Test basic initialization."""
        validator = Validate(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            scoring_functions=[],
        )

        assert str(validator.wrk_dir) == str(temp_dir)
        assert str(validator.system_path) == str(sample_system_yaml)
        assert str(validator.options_path) == str(sample_options_yaml)
        assert validator.sys is not None
        assert validator.opt is not None
        assert validator.scoring_functions == set()
