import os
import logging

def initiate_logger(logger: logging.Logger, debug: bool = False, wrk_dir: str = None, log_name: str = "boltz-eval") -> logging.Logger:
    """
    Set up logging with both file and console handlers.
    
    Args:
        logger (logging.Logger): The logger instance to configure.
        debug (bool, optional): If True, set file logging to DEBUG. Defaults to False.
        wrk_dir (str, optional): Directory to save log file. Defaults to current working directory.
        log_name (str, optional): Base name for log file. Defaults to "boltz-eval".
    
    Returns:
        logging.Logger: Configured logger.
    """
    if wrk_dir is None:
        wrk_dir = os.getcwd()

    log_file = os.path.join(wrk_dir, f"{log_name}.log")
    # Ensure the log file exists and is empty
    open(log_file, 'w+').close()

    # File handler (debug or info)
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG if debug else logging.INFO)

    # Console handler (warnings and above)
    ch = logging.StreamHandler()
    ch.setLevel(logging.WARNING)

    # Formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d,%H:%M:%S'
    )
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    # Clear previous handlers to avoid duplicates
    if logger.hasHandlers():
        logger.handlers.clear()

    logger.setLevel(logging.DEBUG if debug else logging.INFO)
    logger.addHandler(fh)
    logger.addHandler(ch)

    logger.debug(f"Logger initialized. Log file: {log_file}, debug={debug}")
    return logger