#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Docstring____________________________________________________________________

"""**The module *housekeeping_helper* provides general helper functions used
in the ProtBindScreen package.**

Notes:
    The module housekeeping_helper includes functions for:

    - Setting up logging

Imports:
    - Standard libraries: datetime, logging, sys, pathlib

"""

# Imports______________________________________________________________________


import logging
from pathlib import Path
import sys


# Logging utilities___________________________________________________________


def setup_central_logging(
    out_log_path: Path, err_log_path: Path, level: int = logging.INFO
) -> None:
    """Set up the root logger with three handlers:

      - A file handler writing all logs to the `.out` file.
      - A stream handler writing logs to stdout (console).
      - A file handler writing only ERROR and CRITICAL logs to the `.err` file.

    Args:
        out_log_path (Path): Full path to the .out log file.
        err_log_path (Path): Full path to the .err log file.
        level (int): Logging level.

    Side Effects:
        - Removes any pre-existing handlers from the root logger.
        - Installs new handlers on the root logger, affecting all loggers
        in this process.
        - Routes Python warnings into the logging system
        (via `logging.captureWarnings`).
        - Replaces `sys.excepthook` so uncaught exceptions are logged
        automatically.

    Notes:
        This function is intended for CLI entry points, not for library code.
        Because it modifies global logging configuration, it should be called
        once per process.

    """
    # Initialize root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(level)

    # Remove existing handlers to avoid duplicate stream logs
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    # Add new handlers
    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(name)s | %(filename)s:%(lineno)d | "
        "%(funcName)s() | %(message)s"
    )

    # .out file handler
    file_handler = logging.FileHandler(out_log_path)
    file_handler.setFormatter(formatter)
    root_logger.addHandler(file_handler)

    # .stdout stream handler (Console/Jupyter output)
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)
    root_logger.addHandler(stream_handler)

    # .err file for ERROR and CRITICAL logs
    err_handler = logging.FileHandler(err_log_path)
    err_handler.setLevel(logging.ERROR)
    err_handler.setFormatter(formatter)
    root_logger.addHandler(err_handler)

    # Capture warnings and uncaught exceptions
    logging.captureWarnings(True)

    def handle_exception(exc_type, exc_value, exc_traceback):
        if not issubclass(exc_type, KeyboardInterrupt):
            logging.error(
                "Uncaught exception",
                exc_info=(exc_type, exc_value, exc_traceback),
            )

    sys.excepthook = handle_exception

