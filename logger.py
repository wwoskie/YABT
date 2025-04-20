import logging
import random
import string
from pathlib import Path

VERBOSITY_LEVELS = {
    0: logging.WARNING,
    1: logging.INFO,
    2: logging.DEBUG,
}


def setup_class_logger(
    class_name: str,
    log_file: Path | str | None,
    verbosity_console: int,
    verbosiy_file: int,
    logger_format: logging.Formatter = logging.Formatter(
        "%(levelname)s [%(asctime)s] [%(name)s] %(message)s"
    ),
    log_to_console: bool = True,
) -> logging.Logger:
    """Setup unique logger for a class instance"""
    instance_id = "".join(random.choices(string.ascii_uppercase + string.digits, k=6))
    logger_name = f"{class_name}-{instance_id}"
    logger = logging.getLogger(logger_name)

    # if logger.hasHandlers():
    #     return logger

    logger.setLevel(logging.DEBUG)

    if log_to_console:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(VERBOSITY_LEVELS.get(verbosity_console, logging.INFO))
        console_handler.setFormatter(logger_format)
        logger.addHandler(console_handler)

    if log_file is not None:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(VERBOSITY_LEVELS.get(verbosiy_file, logging.INFO))
        file_handler.setFormatter(logger_format)
        logger.addHandler(file_handler)

    return logger
