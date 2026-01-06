"""Tests for boltz_lab.recipes.predict module."""
from unittest.mock import patch


from boltz_lab.recipes.predict import Predict


class TestPredictInit:
    """Tests for Predict initialization."""

    def test_init_basic(self, sample_system_yaml, sample_options_yaml, temp_dir):
        """Test basic initialization."""
        predictor = Predict(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml)
        )

        assert predictor.wrk_dir == str(temp_dir)
        assert predictor.system_path == str(sample_system_yaml)
        assert predictor.options_path == str(sample_options_yaml)
        assert predictor.sys is not None
        assert predictor.opt is not None

    def test_init_with_debug(self, sample_system_yaml, sample_options_yaml, temp_dir):
        """Test initialization with debug logging."""
        predictor = Predict(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            debug=True
        )

        assert predictor.logger.level == 10  # DEBUG level


class TestPredictRun:
    """Tests for Predict.run method."""

    @patch('boltz_lab.recipes.predict.subprocess.run')
    def test_run_executes_command(self, mock_run, sample_system_yaml, sample_options_yaml, temp_dir):
        """Test that run method executes Boltz command."""
        predictor = Predict(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml)
        )

        predictor.run()

        # Verify subprocess.run was called
        assert mock_run.called
        call_args = mock_run.call_args[0][0]
        assert "boltz" in call_args
        assert "predict" in call_args

    @patch('boltz_lab.recipes.predict.subprocess.run')
    def test_run_sets_output_directory(self, mock_run, sample_system_yaml, sample_options_yaml, temp_dir):
        """Test that run method sets output directory."""
        predictor = Predict(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml)
        )

        predictor.run()

        assert predictor.opt.out_dir == str(temp_dir)
