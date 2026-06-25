import marimo

__generated_with = "0.23.8"
app = marimo.App()


@app.cell
def _():
    import sys
    from pathlib import Path

    notebook_dir = Path(__file__).resolve().parent
    if str(notebook_dir) not in sys.path:
        sys.path.insert(0, str(notebook_dir))

    import marimo as mo
    import pandas as pd
    from rdkit import Chem

    from _marimo_helpers import EXAMPLES_DIR, code_block, format_command, make_workspace, read_text

    return Chem, EXAMPLES_DIR, code_block, format_command, make_workspace, mo, pd, read_text


@app.cell
def _(EXAMPLES_DIR, code_block, format_command, mo, read_text):
    covalent_path = EXAMPLES_DIR / "system_covalent.yaml"
    preview_command = [
        "cofolder",
        "validate",
        "-s",
        EXAMPLES_DIR / "system.yaml",
        "-o",
        EXAMPLES_DIR / "options.yaml",
        "--conformers",
        "3D",
        "--runner",
        "boltz2",
    ]
    mo.md(
        f"""
        # Ligand Handling

        This notebook focuses on ligand-specific preparation patterns:

        - choosing ligand input formats
        - understanding COFOLDER conformer flags
        - generating local SDF assets for inspection
        - understanding covalent-system inputs

        **Covalent example system**
        {code_block(read_text(covalent_path), "yaml")}

        **Conformer-enabled command preview**
        ```bash
        {format_command(preview_command)}
        ```

        Notes:

        - `--conformers 2D` produces layout-oriented coordinates
        - `--conformers 3D` produces ETKDG + UFF conformers
        - `--conformers sdf` reuses conformers from an external SDF file
        """
    )
    return


@app.cell
def _(Chem):
    smiles = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
    mol = Chem.MolFromSmiles(smiles)
    canonical_smiles = None if mol is None else Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    return canonical_smiles, smiles


@app.cell
def _(canonical_smiles, mo, smiles):
    mo.md(
        f"""
        ## Quick SMILES Sanity Check

        Example input SMILES:

        ```text
        {smiles}
        ```

        Canonicalized by RDKit:

        ```text
        {canonical_smiles}
        ```
        """
    )
    return


@app.cell
def _(mo):
    run_ligand_demo = mo.ui.checkbox(
        value=False,
        label="Run the local ligand-utility demo (CSV -> SDF -> 2D/3D conformers)",
    )
    run_ligand_demo
    return (run_ligand_demo,)


@app.cell
def _(make_workspace, mo, pd, run_ligand_demo):
    mo.stop(
        not run_ligand_demo.value,
        mo.md(
            "Enable the checkbox above to run the local ligand-utility demo. "
            "This demo assumes a backend-capable environment where "
            "`cofolder.modules.entities.ligand` imports cleanly."
        ),
    )

    try:
        from cofolder.modules.entities.ligand import (
            csv_to_sdf,
            generate_2d_conformers,
            generate_3d_conformers,
            iterate_sdf_records,
            sanitize_mol_id,
        )
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "Ligand utility imports failed. Install a runner environment such as "
            '`python -m pip install -e ".[boltz2]"` before running this demo.'
        ) from exc

    workspace = make_workspace("cofolder-tutorial-ligand-")
    csv_path = workspace / "ligands.csv"
    csv_path.write_text(
        "compound_id,smiles\n"
        "IBUPROFEN,CC(C)Cc1ccc(cc1)C(C)C(=O)O\n"
        "ETHANOL,CCO\n",
        encoding="utf-8",
    )
    sdf_path = workspace / "ligands.sdf"
    sdf_2d_path = workspace / "ligands_2d.sdf"
    sdf_3d_path = workspace / "ligands_3d.sdf"

    csv_to_sdf(
        csv_path=str(csv_path),
        smiles_col="smiles",
        output_sdf_path=str(sdf_path),
        property_cols=["compound_id"],
    )
    generate_2d_conformers(str(sdf_path), str(sdf_2d_path))
    generate_3d_conformers(str(sdf_path), str(sdf_3d_path))
    records = list(iterate_sdf_records(str(sdf_3d_path), "compound_id"))
    summary = pd.DataFrame(
        {
            "compound_id": [record[1] for record in records],
            "sanitized_ccd_candidate": [sanitize_mol_id(record[1]) for record in records],
        }
    )
    return records, sdf_2d_path, sdf_3d_path, summary, workspace


@app.cell
def _(mo, records, sdf_2d_path, sdf_3d_path, summary, workspace):
    mo.md(
        f"""
        ## Utility Demo Results

        Workspace: `{workspace}`

        Generated files:

        - `{sdf_2d_path.name}`
        - `{sdf_3d_path.name}`

        Molecules found in the 3D SDF: `{len(records)}`
        """
    )
    summary
    return


if __name__ == "__main__":
    app.run()
