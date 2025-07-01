# Script containing boltz functions, such as the compiling of the boltz predict
# command.

import multiprocessing

def set_command(yaml_path, i, row):
    command = [
        "boltz",
        "predict",
        yaml_path
    ]

    if i == 1:
        command.extend([f"--use_msa_server"])

    options = opt_data.get("options")
    for n, item in enumerate(options):
        for key, value in item.items():
            value = str(value)

            if key == "out_dir":
                command.extend([f"--out_dir"])
                command.extend([os.path.join(run_dir, f'{i}_{row[id_col]}')])   
                continue     
            
            if key == "use_msa_server":
                continue

            if value == "None" or value == "False":
                continue

            if value == "True":
                command.extend([f"--{key}"])
                continue

            if value == "multiprocessing.cpu_count()":
                command.extend([f"--{key}"])
                command.extend([str(multiprocessing.cpu_count())])
                continue
            
            command.extend([f"--{key}"])
            command.extend([str(value)]) 

    return command