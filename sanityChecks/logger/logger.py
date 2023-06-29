"""
This function will create a new logging file and return a logger object.
This is done because CobraPy can be noisy and does not let us hide information.
"""

import os
import logging
from sanityChecks.utils.utils import get_project_root


def get_logger() -> logging.Logger:
    # Define the logging pattern as "time - level - message"
    # Do not include milliseconds
    logging.basicConfig(
        filename=os.path.join(
            get_project_root(),
            "logs",
            "main.log"
        ),
        level=logging.ERROR,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    logging.getLogger().addHandler(logging.NullHandler())
    return logging.getLogger(__name__)
