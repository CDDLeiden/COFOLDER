import os
import pandas as pd
import subprocess
import time
import json
import shutil
import logging

from rdkit import Chem

from ..helpers import command, conformers, system, utils

def add_arguments(parser):
    """Add screen-specific CLI arguments."""
    parser.add_argument('-w', '--wrk_dir',
                        type=str,
                        dest='wrk_dir',
                        help='Set Working directory if different from CWD.',
                        default=os.getcwd())

    parser.add_argument('-s', '--system_path',
                        type=str,
                        dest='system_path',
                        help='Path to system YAML file.',
                        required=True)

    parser.add_argument('-b', '--boltz_options_path',
                        type=str,
                        dest='options_path',
                        help='Path to Boltz options YAML file.',
                        required=True)
    
    parser.add_argument('-v', '--variable',
                        type=str,
                        dest='variable',
                        help='Comma-seperated list of keys specifying the nested path in \
                            the system YAML to update. (e.g. "sequences,1,ligand,smiles")',
                        required=True)
    
    parser.add_argument('-c', '--variable_csv',
                        type=str,
                        default=None,
                        dest='variable_csv',
                        help='Path to CSV file containing variables. If provided, you must \
                            also specify --col_variable and --col_id.')
    
    parser.add_argument('--col_variable',
                        type=str,
                        default=None,
                        dest='col_variable',
                        help='Column containing variable (e.g. SMILES/CCD/FASTA) (required \
                            if --csv is used).')
    
    parser.add_argument('--col_id',
                        type=str,
                        default=None,
                        dest='col_id',
                        help='Column containing variable ID (required if --csv is used).')
    
    parser.add_argument('-s,', '--variable_sdf',
                        type=str,
                        default=None,
                        dest='variable_sdf',
                        help='Path to SDF file containing variables. If provided, you must \
                            also specify --propterty_id.')
    
    parser.add_argument('--property_id',
                        type=str,
                        default=None,
                        dest='property_id',
                        help='Property name for compound ID in SDF file. Required if \
                            --variable_sdf is used.')
    
    parser.add_argument('--generate_conformers',
                        choices=['2D', '3D'],
                        default=None,
                        dest='generate_conformers',
                        help='Generate 2D or 3D conformers for CCD input. If not specified, \
                            original SMILES (csv) or MolBlock (sdf) are used as system input. \
                            Note: This only works for SMILES, not any other variable type.')
    
    parser.add_argument('--merge_data',
                        type=str,
                        default=None,
                        help='Comma-separated list of CSV columns or SDF prperties from input \
                            to merge into output. (optional)')
    
    parser.add_argument('-d', '--debug',
                        action='store_true',
                        help='Enable debug logging')

def main(args):
    """Run the virtual screening tool."""
    utils.set_dir(args.wrk_dir)  

    logger = logging.getLogger('boltz-tools')
    initiate_logger(logger, debug=args.debug, wrk_dir=args.wrk_dir)
    logger.info("Boltz-tools screen started.")

    screen = Screen(
        wrk_dir=args.wrk_dir,
        system_path=args.system_path,
        options_path=args.options_path,
        variable=args.variable,
        variable_csv=args.variable_csv,
        col_variable=args.col_variable,
        col_id=args.col_id, 
        variable_sdf=args.variable_sdf,
        property_id=args.property_id,
        generate_conformers=args.generate_conformers,
        merge_data=args.merge_data,
    )
    screen.iterate()

def initiate_logger(logger, debug, wrk_dir):
    log_file = os.path.join(wrk_dir, 'boltz-tools.log')
    open(log_file, 'w+').close()
    
    fh = logging.FileHandler(log_file)
    ch = logging.StreamHandler()
    
    fh.setLevel(logging.DEBUG if debug else logging.INFO)
    ch.setLevel(logging.WARNING)
    
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d,%H:%M:%S')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    logger.setLevel(logging.DEBUG if debug else logging.INFO)
    logger.addHandler(fh)
    logger.addHandler(ch)

class Screen(object):
    def __init__(self, **kwargs):
        self.wrk_dir = kwargs.get("wrk_dir")
        self.system_path = kwargs.get("system_path")
        self.options_path = kwargs.get("options_path")
        self._variable = kwargs.get("variable", [])
        self.variable = self.variable = self.variable = [int(v.strip()) if v.strip().isdigit() else v.strip() for v in self._variable.split(",")] if self._variable else []
        
        self.variable_csv = kwargs.get("variable_csv")
        self.col_variable = kwargs.get("col_variable")
        self.col_id = kwargs.get("col_id")
        
        self.variable_sdf = kwargs.get("variable_sdf")
        self.property_id = kwargs.get("property_id")
        
        self.generate_conformers = kwargs.get("generate_conformers")

        self._merge_data = kwargs.get("merge_data", [])
        self.merge_data = [col.strip() for col in self._merge_data.split(",")] if self._merge_data else []
        
        self.logger = logging.getLogger('boltz-tools.screening.Screen')
        self.logger.debug(f"Screen args: {kwargs}")

        # Set command and system objects
        self._options = utils.read_yaml(path=self.options_path)
        self.opt = command.Command(options=self._options)
        self._system = utils.read_yaml(path=self.system_path)
        self.sys = system.System(system=self._system)

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

            # Cache SDF as PKL
            conformers.cache_mols_from_sdf(self.variable_sdf, 
                                            property_id=self.property_id, 
                                            cache=self.opt.find_value(key='cache') or '~/.boltz/')
                
    def iterate(self):      
        if self.variable_sdf:
            for i, mol in enumerate(self.variables, 1):
                start_time = time.time()
                    
                name = mol.GetProp(self.property_id)
                basename = f"{i}_{name}"
                variable = name 

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
