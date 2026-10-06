import logging

from cofolder.cli import main
from cofolder.modules.utils.log import setup_root_logger

def run():
    setup_root_logger(logging.INFO)
    raise SystemExit(main())

if __name__ == "__main__":
    run()
