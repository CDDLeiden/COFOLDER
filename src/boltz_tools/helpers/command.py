# Script containing boltz functions, such as the compiling of the boltz predict
# command.

import multiprocessing
import os

def set_command(yaml_path, opt_data, i, run_dir):
    cmd = [
        "boltz",
        "predict",
        yaml_path
    ]

    cmd.extend([f"--out_dir", run_dir])

    if i == 1:
        cmd.extend([f"--use_msa_server"])

    options = opt_data.get("options")
    for n, item in enumerate(options):
        for key, value in item.items():
            value = str(value)
            
            if key == "use_msa_server":
                continue

            if value == "None" or value == "False":
                continue

            if value == "True":
                cmd.extend([f"--{key}"])
                continue

            if value == "multiprocessing.cpu_count()":
                cmd.extend([f"--{key}"])
                cmd.extend([str(multiprocessing.cpu_count())])
                continue
            
            cmd.extend([f"--{key}"])
            cmd.extend([str(value)]) 

    return cmd