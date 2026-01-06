"""Tests for boltz_lab.modules.utils.helpers module."""

import pytest
import pandas as pd
from rdkit import Chem

from boltz_lab.modules.utils import helpers


class TestCreateDir:
    """Tests for create_dir function."""

    def test_create_new_directory(self, temp_dir):
        """Test creating a new directory."""
        new_dir = temp_dir / "new_test_dir"
        assert not new_dir.exists()

        helpers.create_dir(str(new_dir))

        assert new_dir.exists()
        assert new_dir.is_dir()

    def test_existing_directory(self, temp_dir):
        """Test with an existing directory."""
        helpers.create_dir(str(temp_dir))
        assert temp_dir.exists()


class TestSetDir:
    """Tests for set_dir function."""

    def test_create_new_directory(self, temp_dir):
        """Test creating a new directory."""
        new_dir = temp_dir / "set_test_dir"
        assert not new_dir.exists()

        helpers.set_dir(str(new_dir))

        assert new_dir.exists()
        assert new_dir.is_dir()


class TestParseListAsStr:
    """Tests for parse_list_as_str function."""

    def test_parse_integers(self):
        """Test parsing comma-separated integers."""
        result = helpers.parse_list_as_str("1,2,3", item_type=int)
        assert result == [1, 2, 3]

    def test_parse_strings(self):
        """Test parsing comma-separated strings."""
        result = helpers.parse_list_as_str("a,b,c", item_type=str)
        assert result == ["a", "b", "c"]

    def test_parse_floats(self):
        """Test parsing comma-separated floats."""
        result = helpers.parse_list_as_str("1.5,2.3,3.7", item_type=float)
        assert result == [1.5, 2.3, 3.7]

    def test_custom_separator(self):
        """Test parsing with custom separator."""
        result = helpers.parse_list_as_str("1;2;3", separator=";", item_type=int)
        assert result == [1, 2, 3]

    def test_expected_length(self):
        """Test with expected length validation."""
        result = helpers.parse_list_as_str("1,2,3", item_type=int, expected_length=3)
        assert result == [1, 2, 3]

    def test_wrong_expected_length(self):
        """Test with wrong expected length."""
        with pytest.raises(ValueError, match="List length must be"):
            helpers.parse_list_as_str("1,2,3", item_type=int, expected_length=2)

    def test_invalid_type_conversion(self):
        """Test with invalid type conversion."""
        with pytest.raises(ValueError, match="Failed to convert"):
            helpers.parse_list_as_str("a,b,c", item_type=int)

    def test_strip_whitespace(self):
        """Test that whitespace is stripped."""
        result = helpers.parse_list_as_str(" 1 , 2 , 3 ", item_type=int)
        assert result == [1, 2, 3]


class TestReadYaml:
    """Tests for read_yaml function."""

    def test_read_valid_yaml(self, sample_yaml_file):
        """Test reading a valid YAML file."""
        data = helpers.read_yaml(str(sample_yaml_file))

        assert data is not None
        assert data["key1"] == "value1"
        assert data["key2"]["nested"] == "value2"

    def test_read_nonexistent_file(self, temp_dir):
        """Test reading a nonexistent file."""
        result = helpers.read_yaml(str(temp_dir / "nonexistent.yaml"))
        assert result is None


class TestReadCsv:
    """Tests for read_csv function."""

    def test_read_valid_csv(self, sample_csv_file):
        """Test reading a valid CSV file."""
        df = helpers.read_csv(str(sample_csv_file), columns=["compound_id", "smiles"])

        assert df is not None
        assert len(df) == 2
        assert "compound_id" in df.columns
        assert "smiles" in df.columns

    def test_missing_column_warning(self, sample_csv_file, caplog):
        """Test warning for missing column."""
        df = helpers.read_csv(str(sample_csv_file), columns=["nonexistent"])
        assert df is not None


class TestReadSdf:
    """Tests for read_sdf function."""

    def test_read_valid_sdf(self, temp_dir):
        """Test reading a valid SDF file."""
        sdf_path = temp_dir / "test.sdf"

        # Create a simple SDF file
        writer = Chem.SDWriter(str(sdf_path))
        mol1 = Chem.MolFromSmiles("CCO")
        mol2 = Chem.MolFromSmiles("CC(C)O")
        writer.write(mol1)
        writer.write(mol2)
        writer.close()

        mols = helpers.read_sdf(str(sdf_path))

        assert mols is not None
        assert len(mols) == 2


class TestDeleteLastLine:
    """Tests for delete_last_line function."""

    def test_delete_last_line(self, temp_dir):
        """Test deleting the last line of a file."""
        test_file = temp_dir / "test.txt"
        test_file.write_text("line1\nline2\nline3\n")

        helpers.delete_last_line(str(test_file))

        content = test_file.read_text()
        assert content == "line1\nline2\n"

    def test_delete_from_single_line(self, temp_dir):
        """Test deleting from a single-line file."""
        test_file = temp_dir / "test.txt"
        test_file.write_text("single line\n")

        helpers.delete_last_line(str(test_file))

        content = test_file.read_text()
        assert content == ""

    def test_empty_file(self, temp_dir):
        """Test with an empty file."""
        test_file = temp_dir / "test.txt"
        test_file.write_text("")

        helpers.delete_last_line(str(test_file))

        content = test_file.read_text()
        assert content == ""


class TestParseCensoredAffinity:
    """Tests for parse_censored_affinity function."""

    def test_parse_with_signs(self):
        """Test parsing affinity values with censoring signs."""
        series = pd.Series([">5.0", "<=6.5", "7.2", "<4.0"])
        result = helpers.parse_censored_affinity(series, keep_sign=True)

        assert result["affinity_value"].tolist() == [5.0, 6.5, 7.2, 4.0]
        assert result["affinity_sign"].tolist() == [">", "<=", None, "<"]

    def test_parse_without_signs(self):
        """Test parsing with signs removed."""
        series = pd.Series([">5.0", "6.5"])
        result = helpers.parse_censored_affinity(series, keep_sign=False)

        assert all(pd.isnull(result["affinity_sign"]))


class TestRemoveCensoredAffinity:
    """Tests for remove_censored_affinity function."""

    def test_remove_censored_rows(self):
        """Test removing rows with censoring signs."""
        df = pd.DataFrame({
            "affinity": [">5.0", "6.5", "<=7.2", "8.0"],
            "other": [1, 2, 3, 4]
        })
        result = helpers.remove_censored_affinity(df, ["affinity"])

        assert len(result) == 2
        assert result["affinity"].tolist() == ["6.5", "8.0"]


class TestStripCensoringSigns:
    """Tests for strip_censoring_signs function."""

    def test_strip_signs(self):
        """Test stripping censoring signs."""
        df = pd.DataFrame({
            "aff1": [">5.0", "6.5", "<=7.2"],
            "aff2": ["<4.0", "5.5", ">=6.0"]
        })
        result = helpers.strip_censoring_signs(df, ["aff1", "aff2"])

        assert result["aff1"].tolist() == [5.0, 6.5, 7.2]
        assert result["aff2"].tolist() == [4.0, 5.5, 6.0]


class TestPrepareAffinityDataframe:
    """Tests for prepare_affinity_dataframe function."""

    def test_remove_mode(self):
        """Test with censoring='remove'."""
        df = pd.DataFrame({
            "pred": [">5.0", "6.5", "7.0"],
            "exp": ["5.5", "6.0", "<7.5"]
        })
        result = helpers.prepare_affinity_dataframe(df, ["pred", "exp"], censoring="remove")

        assert len(result) == 1
        assert result["pred"].tolist() == ["6.5"]

    def test_strip_mode(self):
        """Test with censoring='strip'."""
        df = pd.DataFrame({
            "pred": [">5.0", "6.5"],
            "exp": ["5.5", "<6.0"]
        })
        result = helpers.prepare_affinity_dataframe(df, ["pred", "exp"], censoring="strip")

        assert len(result) == 2
        assert result["pred"].tolist() == [5.0, 6.5]
        assert result["exp"].tolist() == [5.5, 6.0]

    def test_invalid_mode(self):
        """Test with invalid censoring mode."""
        df = pd.DataFrame({"pred": [5.0], "exp": [5.0]})
        with pytest.raises(ValueError, match="censoring must be"):
            helpers.prepare_affinity_dataframe(df, ["pred", "exp"], censoring="invalid")


class TestConvertBoltzAffinityToIc50:
    """Tests for convert_boltz_affinity_to_ic50 function."""

    def test_conversion(self):
        """Test affinity conversion."""
        df = pd.DataFrame({
            "affinity_pred_value": [6.0, 7.0, 8.0]
        })
        result = helpers.convert_boltz_affinity_to_ic50(df)

        assert "IC50_uM" in result.columns
        assert "pIC50_kcal_per_mol" in result.columns
        assert result["IC50_uM"].tolist() == [1e6, 1e7, 1e8]

    def test_with_output_path(self, temp_dir):
        """Test with output file."""
        df = pd.DataFrame({"affinity_pred_value": [6.0]})
        output_path = temp_dir / "output.csv"

        helpers.convert_boltz_affinity_to_ic50(df, output_path=str(output_path))

        assert output_path.exists()


class TestDropAndLogNans:
    """Tests for drop_and_log_nans function."""

    def test_drop_nans(self):
        """Test dropping rows with NaN values."""
        df = pd.DataFrame({
            "a": [1.0, 2.0, float("nan"), 4.0],
            "b": [5.0, float("nan"), 7.0, 8.0],
            "c": [9.0, 10.0, 11.0, 12.0]
        })
        result = helpers.drop_and_log_nans(df, ["a", "b"])

        assert len(result) == 2
        assert result["a"].tolist() == [1.0, 4.0]


class TestCalculateAffinityCorrelations:
    """Tests for calculate_affinity_correlations function."""

    def test_calculate_correlations(self):
        """Test correlation calculation."""
        df = pd.DataFrame({
            "pred": [5.0, 6.0, 7.0, 8.0, 9.0],
            "exp": [5.2, 6.1, 6.9, 8.2, 8.8]
        })
        metrics = helpers.calculate_affinity_correlations(df, "pred", "exp")

        assert "r2" in metrics
        assert "pearson" in metrics
        assert "spearman" in metrics
        assert "kendall" in metrics
        assert "rmse" in metrics
        assert "mae" in metrics

        # Check that correlations are reasonable
        assert 0 <= metrics["r2"] <= 1
        assert -1 <= metrics["pearson"] <= 1
        assert -1 <= metrics["spearman"] <= 1


class TestPlotAffinityCorrelation:
    """Tests for plot_affinity_correlation function."""

    def test_plot_to_file(self, temp_dir):
        """Test plotting to a file."""
        df = pd.DataFrame({
            "pred": [5.0, 6.0, 7.0, 8.0],
            "exp": [5.2, 6.1, 6.9, 8.2]
        })
        output_path = temp_dir / "plot.png"

        helpers.plot_affinity_correlation(
            df, "pred", "exp", output_path=str(output_path)
        )

        assert output_path.exists()
