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
        assert validator.reference_path is None
        assert validator.reproduction_metrics == {
            "protein_rmsd",
            "ligand_rmsd",
            "sucos",
            "pocket_coverage",
        }

    def test_init_with_reference_path(self, sample_system_yaml, sample_options_yaml, temp_dir):
        """Test initialization with a valid reference structure path."""
        reference_path = temp_dir / "reference.pdb"
        reference_path.write_text("HEADER TEST\n")

        validator = Validate(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            scoring_functions=[],
            reference_path=str(reference_path),
            reproduction_metrics=["sucos"],
        )

        assert str(validator.reference_path) == str(reference_path)
        assert validator.reproduction_metrics == {"sucos"}
