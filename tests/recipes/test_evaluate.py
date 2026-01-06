"""Tests for boltz_lab.recipes.evaluate module."""
from unittest.mock import Mock

import pytest

from boltz_lab.recipes.evaluate import Evaluate


class TestEvaluateInit:
    """Tests for Evaluate initialization."""

    def test_init_basic(self, sample_system_yaml, sample_options_yaml, temp_dir):
        """Test basic initialization."""
        evaluator = Evaluate(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml)
        )

        assert evaluator.wrk_dir == str(temp_dir)
        assert evaluator.system_path == str(sample_system_yaml)
        assert evaluator.options_path == str(sample_options_yaml)
        assert evaluator.repeats == 1
        assert evaluator.seeds == []

    def test_init_with_repeats(self, sample_system_yaml, sample_options_yaml, temp_dir):
        """Test initialization with repeats."""
        evaluator = Evaluate(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            repeats=5
        )

        assert evaluator.repeats == 5

    def test_init_with_seeds(self, sample_system_yaml, sample_options_yaml, temp_dir):
        """Test initialization with seeds."""
        evaluator = Evaluate(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            repeats=3,
            seeds="1,2,3"
        )

        assert evaluator.seeds == [1, 2, 3]

    def test_init_with_input_pdb(self, sample_system_yaml, sample_options_yaml, temp_dir):
        """Test initialization with input PDB."""
        evaluator = Evaluate(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            input_pdb="/path/to/structure.pdb"
        )

        assert evaluator.input_pdb == "/path/to/structure.pdb"


class TestEvaluateLoadReferenceStructure:
    """Tests for Evaluate._load_reference_structure method."""

    def test_load_reference_none(self, sample_system_yaml, sample_options_yaml, temp_dir):
        """Test loading reference structure when none provided."""
        evaluator = Evaluate(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml)
        )

        result = evaluator._load_reference_structure()

        assert result is None

    def test_load_reference_unsupported_format(self, sample_system_yaml, sample_options_yaml, temp_dir):
        """Test loading unsupported file format raises error."""
        evaluator = Evaluate(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            input_pdb="/path/to/structure.xyz"
        )

        with pytest.raises(RuntimeError, match="Unsupported file type"):
            evaluator._load_reference_structure()


class TestEvaluateExtractIFP:
    """Tests for Evaluate._extract_ifp method."""

    def test_extract_ifp_false(self, sample_system_yaml, sample_options_yaml, temp_dir):
        """Test IFP extraction when ifp=false."""
        evaluator = Evaluate(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            ifp="false"
        )

        evaluator._extract_ifp(Mock())

        assert evaluator.ifp_data is None

    def test_extract_ifp_none(self, sample_system_yaml, sample_options_yaml, temp_dir):
        """Test IFP extraction when ifp is None."""
        evaluator = Evaluate(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml),
            ifp=None
        )

        evaluator._extract_ifp(Mock())

        assert evaluator.ifp_data is None


class TestEvaluateParseSeeds:
    """Tests for Evaluate._parse_seeds method."""

    def test_parse_seeds_valid(self, sample_system_yaml, sample_options_yaml, temp_dir):
        """Test parsing valid seeds string."""
        evaluator = Evaluate(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml)
        )

        result = evaluator._parse_seeds("1,2,3")

        assert result == [1, 2, 3]

    def test_parse_seeds_none(self, sample_system_yaml, sample_options_yaml, temp_dir):
        """Test parsing None seeds."""
        evaluator = Evaluate(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml)
        )

        result = evaluator._parse_seeds(None)

        assert result == []

    def test_parse_seeds_empty(self, sample_system_yaml, sample_options_yaml, temp_dir):
        """Test parsing empty seeds string."""
        evaluator = Evaluate(
            wrk_dir=str(temp_dir),
            system_path=str(sample_system_yaml),
            options_path=str(sample_options_yaml)
        )

        result = evaluator._parse_seeds("")

        assert result == []
