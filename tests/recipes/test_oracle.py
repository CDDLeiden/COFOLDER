"""Tests for cofolder.recipes.oracle module."""
import logging


from cofolder.recipes.oracle import Oracle


class TestOracleInit:
    """Tests for Oracle initialization."""

    def test_init_basic(self, sample_system_yaml, sample_options_yaml, temp_dir):
        """Test basic initialization."""
        oracle = Oracle(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml)
        )

        assert oracle.wrk_dir == str(temp_dir)
        assert oracle.system_path == str(sample_system_yaml)
        assert oracle.options_path == str(sample_options_yaml)

    def test_init_with_debug(self, sample_system_yaml, sample_options_yaml, temp_dir):
        """Test initialization with debug enabled."""
        oracle = Oracle(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            debug=True
        )

        assert oracle.logger.level == logging.DEBUG

    def test_init_without_debug(self, sample_system_yaml, sample_options_yaml, temp_dir):
        """Test initialization without debug (default)."""
        oracle = Oracle(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            debug=False
        )

        assert oracle.logger.level == logging.INFO
