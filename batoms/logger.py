"""
Logger
"""
import bpy
import sys
from tempfile import gettempdir
from pathlib import Path
import logging

logger = logging.getLogger('batoms')

def set_logger():
    formatter = ('%(levelname)-8s '
                                '[%(funcName)-20s]: %(message)s')
    logging.basicConfig(stream=sys.stdout,
                        format=formatter,
                        level=logging.INFO
                        )
    # add logger file
    filepath = Path(gettempdir()) / ("beautiful_atoms.log")
    logger.info("Log file: " + str(filepath))
    file_handler = logging.FileHandler(filepath, mode="w")
    file_handler.setFormatter(logging.Formatter(formatter))
    logger.addHandler(file_handler)

def update_logging_level():
    prefs = bpy.context.preferences.addons['batoms'].preferences
    logger.setLevel(prefs.logging_level)


def print_time(key, value):
    return '{}: {}'.format(key, value)