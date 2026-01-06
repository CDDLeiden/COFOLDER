import logging

from boltz_lab.cli import main
from boltz_lab.modules.utils.log import setup_root_logger

def run():
    setup_root_logger(logging.INFO)
    main()

if __name__ == "__main__":
    run()
    