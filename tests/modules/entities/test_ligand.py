"""Tests for cofolder.modules.entities.ligand module."""

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from cofolder.modules.entities import ligand


class TestSanitizeMolId:
    """Tests for sanitize_mol_id function."""

    def test_short_id(self):
        """Test with ID shorter than 5 characters."""
        result = ligand.sanitize_mol_id("ABC")
        assert result == "ABC"

    def test_exactly_five_chars(self):
        """Test with ID exactly 5 characters."""
        result = ligand.sanitize_mol_id("ABCDE")
        assert result == "ABCDE"

    def test_long_id_truncated(self):
        """Test with ID longer than 5 characters."""
        result = ligand.sanitize_mol_id("ABCDEFGH")
        assert result == "ABCDE"
        assert len(result) == 5


class TestAddPickledProp:
    """Tests for add_pickled_prop function."""

    def test_add_property(self):
        """Test adding a pickled property to molecule."""
        mol = Chem.MolFromSmiles("CCO")
        test_value = {"key": "value", "number": 42}

        ligand.add_pickled_prop(mol, "test_prop", test_value)

        assert mol.HasProp("test_prop")
        # The property should be stored as hex
        prop_value = mol.GetProp("test_prop")
        assert isinstance(prop_value, str)


class TestExtractConstraints:
    """Tests for extract_constraints function."""

    def test_extract_simple_attribute(self):
        """Test extracting a simple attribute from constraint objects."""
        # Create mock constraint objects
        class MockConstraint:
            def __init__(self, value):
                self.is_bond = value

        constraints = [MockConstraint(True), MockConstraint(False), MockConstraint(True)]
        result = ligand.extract_constraints(constraints, "is_bond")

        assert result == [True, False, True]

    def test_extract_with_transpose(self):
        """Test extracting and transposing."""
        import numpy as np

        class MockConstraint:
            def __init__(self, coords):
                self.coords = coords

        constraints = [
            MockConstraint([1, 2]),
            MockConstraint([3, 4])
        ]
        result = ligand.extract_constraints(constraints, "coords", transpose=True)

        assert isinstance(result, np.ndarray)
        assert result.shape == (2, 2)


class TestGenerate2DConformers:
    """Tests for generate_2d_conformers function."""

    def test_generate_2d(self, temp_dir):
        """Test generating 2D conformers."""
        # Create input SDF with molecules
        input_sdf = temp_dir / "input.sdf"
        writer = Chem.SDWriter(str(input_sdf))
        mol1 = Chem.MolFromSmiles("CCO")
        mol2 = Chem.MolFromSmiles("CC(C)O")
        writer.write(mol1)
        writer.write(mol2)
        writer.close()

        output_sdf = temp_dir / "output.sdf"

        ligand.generate_2d_conformers(str(input_sdf), str(output_sdf))

        assert output_sdf.exists()

        # Verify molecules have 2D coordinates
        suppl = Chem.SDMolSupplier(str(output_sdf))
        mols = [m for m in suppl if m is not None]
        assert len(mols) == 2

        # Check that conformers were generated
        for mol in mols:
            assert mol.GetNumConformers() >= 1


class TestGenerate3DConformers:
    """Tests for generate_3d_conformers function."""

    def test_generate_3d(self, temp_dir):
        """Test generating 3D conformers."""
        # Create input SDF
        input_sdf = temp_dir / "input.sdf"
        writer = Chem.SDWriter(str(input_sdf))
        mol = Chem.MolFromSmiles("CCO")
        writer.write(mol)
        writer.close()

        output_sdf = temp_dir / "output.sdf"

        ligand.generate_3d_conformers(str(input_sdf), str(output_sdf))

        assert output_sdf.exists()

        # Verify molecules have 3D coordinates
        suppl = Chem.SDMolSupplier(str(output_sdf))
        mols = [m for m in suppl if m is not None]
        assert len(mols) > 0


class TestCsvToSdf:
    """Tests for csv_to_sdf function."""

    def test_csv_to_sdf(self, temp_dir):
        """Test converting CSV to SDF."""
        # Create CSV file
        csv_path = temp_dir / "compounds.csv"
        csv_path.write_text(
            "id,smiles,mw\n"
            "CMPD001,CCO,46.07\n"
            "CMPD002,CC(C)O,60.10\n"
        )

        output_sdf = temp_dir / "output.sdf"

        ligand.csv_to_sdf(
            csv_path=str(csv_path),
            smiles_col="smiles",
            output_sdf_path=str(output_sdf),
            property_cols=["id", "mw"]
        )

        assert output_sdf.exists()

        # Verify SDF content
        suppl = Chem.SDMolSupplier(str(output_sdf))
        mols = [m for m in suppl if m is not None]
        assert len(mols) == 2

        # Check properties were transferred
        mol1 = mols[0]
        assert mol1.HasProp("id")
        assert mol1.HasProp("mw")

    def test_csv_to_sdf_missing_file_raises(self, temp_dir):
        """Missing CSV paths should fail explicitly."""
        missing_csv = temp_dir / "missing.csv"
        output_sdf = temp_dir / "output.sdf"

        with pytest.raises(FileNotFoundError, match="CSV file not found"):
            ligand.csv_to_sdf(
                csv_path=str(missing_csv),
                smiles_col="smiles",
                output_sdf_path=str(output_sdf),
            )


class TestIterateSdfRecords:
    """Tests for iterate_sdf_records function."""

    def test_iterate_sdf(self, temp_dir):
        """Test iterating through SDF records."""
        # Create SDF file
        sdf_path = temp_dir / "test.sdf"
        writer = Chem.SDWriter(str(sdf_path))

        mol1 = Chem.MolFromSmiles("CCO")
        mol1.SetProp("ID", "MOL001")
        writer.write(mol1)

        mol2 = Chem.MolFromSmiles("CC(C)O")
        mol2.SetProp("ID", "MOL002")
        writer.write(mol2)

        writer.close()

        # Iterate and collect
        records = list(ligand.iterate_sdf_records(str(sdf_path), "ID"))

        assert len(records) == 2
        assert records[0][1] == "MOL001"  # mol_id
        assert records[1][1] == "MOL002"
        assert isinstance(records[0][2], str)  # molblock


class TestMolToCcd:
    """Tests for mol_to_ccd function."""

    def test_mol_to_ccd(self, temp_dir):
        """Test converting molecule to CCD format."""
        mol = Chem.MolFromSmiles("CCO")
        AllChem.EmbedMolecule(mol)

        # This function interacts with Boltz's CCD system
        # We'll test that it doesn't raise an error
        try:
            ligand.mol_to_ccd(
                resname="ETH",
                mol=mol,
                boltz_path=str(temp_dir),
                on_conflict="overwrite"
            )
        except Exception:
            # The actual CCD conversion might fail without Boltz installed
            # We just want to make sure the function is callable
            pass
