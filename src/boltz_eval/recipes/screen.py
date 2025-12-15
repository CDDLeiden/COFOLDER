import os
import pandas as pd
import subprocess
import time
import json
import shutil
import logging

# legacy imports
from ..modules.legacy import command, system, utils

class Screen:
    """High-level orchestrator for virtual screening workflow."""
    def __init__(
        self,
        wrk_dir: str,
        system_path: str,
        options_path: str,
        variable: str | None = None,
        variable_csv: str | None = None,
        col_variable: str | None = None,
        col_id: str | None = None,
        variable_sdf: str | None = None,
        property_id: str | None = None,
        generate_conformers: str | None = None,
        merge_data: str | None = None,
        debug: bool = False,
    ):
        self.wrk_dir = wrk_dir
        self.system_path = system_path
        self.options_path = options_path

        # Parse variables (comma-separated string)
        self.variable = self._parse_list(variable)
        self.merge_data = self._parse_list(merge_data)

        self.variable_csv = variable_csv
        self.col_variable = col_variable
        self.col_id = col_id
        self.variable_sdf = variable_sdf
        self.property_id = property_id
        self.generate_conformers = generate_conformers

        # Logger setup
        self.logger = logging.getLogger('boltz-eval.screening.Screen')
        self.logger.setLevel(logging.DEBUG if debug else logging.INFO)
        self.logger.debug("Initializing Screen with parameters: %s", {
            "wrk_dir": wrk_dir,
            "system_path": system_path,
            "options_path": options_path,
            "variable": variable,
            "variable_csv": variable_csv,
            "col_variable": col_variable,
            "col_id": col_id,
            "variable_sdf": variable_sdf,
            "property_id": property_id,
            "generate_conformers": generate_conformers,
            "merge_data": merge_data,
        })

        # Load system and options
        self._options = utils.read_yaml(path=self.options_path)
        self.opt = command.Command(options=self._options)

        self._system = utils.read_yaml(path=self.system_path)
        self.sys = system.System(system=self._system)

        self.logger.debug("Screen initialization complete.")

    def run(self):
        self.load_screen()
        self.iterate()

    def load_screen(self):
        # Load variables from CSV if provided
        if self.variable_csv:
            if not self.col_variable or not self.col_id:
                raise ValueError("Both col_smiles and col_id must be specified.")

            if self.generate_conformers:
                self.variable_sdf = os.path.splitext(self.variable_csv)[0] + ".sdf"
                self.property_id = self.col_id
                conformers.csv_to_sdf(
                    csv_path=self.variable_csv,
                    smiles_col=self.col_variable,
                    output_sdf_path=self.variable_sdf,
                    property_cols=[self.col_id]+self.merge_data
                )     

            else:
                self.variables = utils.read_csv(
                    path=self.variable_csv, 
                    columns=[self.col_variable, self.col_id] + self.merge_data
                ) 

                missing = [c for c in [self.col_variable, self.col_id] + self.merge_data 
                        if c not in self.variables.columns]
                if missing:
                    raise ValueError(f"CSV missing columns: {', '.join(missing)}")
        
        # If SDF exists (either provided or generated from CSV)
        if self.variable_sdf:
            self.variables = utils.read_sdf(path=self.variable_sdf)
            #self.variables = conformers.refine_sdf(self._variables)

            if self.generate_conformers == "2D":
                conformers.generate_2d_conformers(self.variable_sdf)
            elif self.generate_conformers == "3D":
                conformers.generate_3d_conformers(self.variable_sdf)

            # Cache each mol as PKL
            cache_dir = self.opt.find_value(key='cache') or '~/.boltz/'
            for mol in self.variables:
                if mol is None:
                    continue
                if mol.HasProp(self.property_id):
                    mol_id = mol.GetProp(self.property_id)
                    mol_id = conformers.sanitize_mol_id(mol_id)
                else:
                    self.logger.error(f"SDF molecule missing ID property '{self.property_id}'")
                    continue
                try:
                    conformers.mol_to_ccd(mol_id, mol, boltz_path=cache_dir)
                except Exception as e:
                    self.logger.error(f"Failed to process ID {mol_id}: {e}")
                
    def iterate(self):      
        if self.variable_sdf:
            for i, mol in enumerate(self.variables, 1):
                start_time = time.time()
                    
                name = mol.GetProp(self.property_id)
                basename = f"{i}_{name}"

                variable = name
                if len(name) > 5:
                    truncated = name[:5]
                    print(f"[WARNING] Variable name '{name}' is longer than 5 characters. "
                        f"Truncating to '{truncated}' to comply with CCD naming rules.")
                    variable = truncated


                self.run(variable, basename)

                self.logger.info(" pred time--- %.2f seconds ---" % (time.time() - start_time))

        # Set basename and variable for CSV input         
        elif self.variable_csv:
            # Run the screening process for each variable
            for i, (_, row) in enumerate(self.variables.iterrows(), 1):
                start_time = time.time()

                name = str(row[self.col_id])
                basename = f'{i}_{name}'
                variable = row[self.col_variable]
                self.logger.info(f"({i}/{len(self.variables)}) {name}: {variable}")

                self.run(variable, basename)
                
                self.logger.info(" pred time--- %.2f seconds ---" % (time.time() - start_time))

    def run(self, variable, basename):       
        self.sys.update_system(value=variable, path=self.variable)
        
        # Set output directory and update system
        out_dir = os.path.join(self.wrk_dir, basename)
        utils.set_dir(out_dir)
        self.opt.out_dir = out_dir
        
        # MSA recycling - only possible for monomer systems
        self.msa = self.sys.find_value(key='msa')
        if self.msa and not os.path.exists(self.msa):
            raise FileNotFoundError(f"MSA file not found at {self.msa}")

        self.logger.info(f'System:\n{str(self.sys.system)}')
        yaml_path = os.path.join(out_dir, f'{basename}.yaml')
        self.opt.system_path = yaml_path
        self.sys.save_system_to_yaml(path=yaml_path)

        # Set and run command
        cmd = self.opt.set_command(system=self.sys)
        self.logger.info(f'Running: {" ".join(cmd)}')
        subprocess.run(cmd)

        # update MSA after first iteration
        if self.msa == None:
            try:
                _ = self.sys.find_value(key='protein')
                self.msa = os.path.join(out_dir, f'boltz_results_{basename}/msa/{basename}_unpaired_tmp_env/uniref.a3m')
                self.sys.update_system(value=self.msa, parent_key='protein', sub_key='msa')
                self.logger.info(f'Cleaning up MSA file: {self.msa}')
                utils.delete_last_line(self.msa)
            except ValueError:
                self.logger.info(f'MSA recycling not available for multimers in current version')

        #TODO: gather results (sdf/csv independent) | Current: Only CSV
        #self.gather_metrics(out_dir, i, row)
        #self.gather_structures(out_dir, i, row)


    def gather_metrics(self, out_dir, i, row):
        # function for gathering confidence and affinity metrics into csv
        basename = f'{i}_{row[self.col_id]}'
        csv_file = os.path.join(self.wrk_dir, "output.csv")
        conf_file = os.path.join(out_dir, f'boltz_results_{basename}/predictions/{basename}/confidence_{basename}_model_0.json')
        aff_file = os.path.join(out_dir, f'boltz_results_{basename}/predictions/{basename}/affinity_{basename}.json')

        if not os.path.exists(conf_file):
            self.logger.warning(f"Missing confidence file: {conf_file}")
            return

        try:
            df_conf = pd.json_normalize(json.load(open(conf_file)))
        except Exception as e:
            self.logger.error(f"Error reading {conf_file}: {e}")
            return

        df_new = df_conf
        if os.path.exists(aff_file):
            try:
                df_aff = pd.json_normalize(json.load(open(aff_file)))
                df_new = pd.concat([df_aff, df_conf], axis=1)
            except Exception as e:
                self.logger.error(f"Error reading {aff_file}: {e}")

        output = {**{'index': i, 'id': row[self.col_id], 'basename': basename, 'smiles': row[self.col_variable]},
                  **{col: row[col] for col in self.merge_data},
                  **df_new.iloc[0].to_dict()}

        pd.DataFrame([output]).to_csv(csv_file, mode='a', header=not os.path.exists(csv_file), index=False)

    def gather_structures(self, out_dir, i, row):
        # function to gather cif or pdb files into single folder
        basename = f'{i}_{row[self.col_id]}'
        target_dir = os.path.join(self.wrk_dir, "structures")
        os.makedirs(target_dir, exist_ok=True)
        structure_dir = os.path.join(out_dir, f'boltz_results_{basename}/predictions/{basename}/')
        if not os.path.exists(structure_dir):
            self.logger.warning(f"Missing structure directory: {structure_dir}")
            return
        for file_name in os.listdir(structure_dir):
            if file_name.endswith(('.cif', '.pdb')):
                try:
                    shutil.copy2(os.path.join(structure_dir, file_name), os.path.join(target_dir, file_name))
                except Exception as e:
                    self.logger.error(f"Error copying {file_name}: {e}")

    @staticmethod
    def _parse_list(input_str: str | None) -> list:
        """Parse a comma-separated string into a list, converting digits to int."""
        if not input_str:
            return []
        return [int(v.strip()) if v.strip().isdigit() else v.strip() for v in input_str.split(",")]