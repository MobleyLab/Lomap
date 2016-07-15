r"""
Logger setup (default is no output) and custom logging formatter.
"""

import os
import logging


class LogFormatter(logging.Formatter):
    """
    Custom logging formatter.
    """

    err_fmt  = "ERROR: %(msg)s"
    warn_fmt = "WARNING: %(msg)s"
    info_fmt = "%(msg)s"
    dbg_fmt  = "DEBUG: %(module)s: %(lineno)d: %(msg)s"

    def __init__(self, fmt="%(levelno)s: %(msg)s"):
        logging.Formatter.__init__(self, fmt)

    def format(self, record):
        format_orig = self._fmt

        if record.levelno == logging.DEBUG:
            self._fmt = LogFormatter.dbg_fmt
        elif record.levelno == logging.WARN:
            self._fmt = LogFormatter.warn_fmt
        elif record.levelno == logging.INFO:
            self._fmt = LogFormatter.info_fmt
        elif record.levelno == logging.ERROR:
            self._fmt = LogFormatter.err_fmt

        result = logging.Formatter.format(self, record)

        self._fmt = format_orig

        return result


logging.basicConfig(level=logging.DEBUG, filename=os.devnull)
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

