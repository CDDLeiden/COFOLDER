from boltz_eval.modules.utils.log import setup_root_logger
from boltz_eval.cli import main

import logging

def run():
    setup_root_logger(logging.INFO)
    main()

if __name__ == "__main__":
    run()