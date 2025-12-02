import logging

class Oracle(object):
    """High-level orchestrator for oracle workflow."""
    def __init__(
        self,
        wrk_dir: str,
        system_path: str,
        options_path: str,
        debug: bool = False,
    ):
        self.wrk_dir = wrk_dir
        self.system_path = system_path
        self.options_path = options_path

        # Setup logger
        self.logger = logging.getLogger(f'boltz-eval.oracle.Oracle')
        self.logger.setLevel(logging.DEBUG if debug else logging.INFO)
        self.logger.debug("Initializing Oracle with parameters: %s", {
            "wrk_dir": wrk_dir,
            "system_path": system_path,
            "options_path": options_path
        })

        self.logger.debug("Oracle initialization complete.")