from __future__ import annotations

"""
MCRT2 Python utilities prompt
"""

from collections.abc import Iterable
import sys
import six


def warning(msg: str) -> str:
    return "[WARNING] {}".format(msg)


def error(msg: str) -> str:
    return "[ERROR] {}".format(msg)


def info(msg: str) -> str:
    return "[INFO] {}".format(msg)


class Prompt:
    def __init__(self, dim: int = 1):
        self._dim = dim
    def __getitem__(self, args: Iterable | str):
        n = len(sys.argv)
        if isinstance(args, six.string_types):
            item_list = [args]
        elif isinstance(args, Iterable):
            item_list = args
        else:
            item_list = [args]
        for i in range(1, n):
            arg = sys.argv[i]
            for item in item_list:
                if arg == item:
                    return sys.argv[i + 1: i + 1 + self._dim] if i + self._dim < n else ""
        return None
