from __future__ import annotations

"""
ENDF Flags

this code is part of the RT2 project
"""

from src.prompt import warning
from enum import Enum


class ZSYMAM:
    def __init__(self, zsymam_str: str, verbose: bool):
        self._repr = zsymam_str
        try:
            zstr = zsymam_str[0:3]
            astr = zsymam_str[7:10]

            self._ss   = zsymam_str[4:6]
            self._z    = 0 if zstr.isspace() else int(zstr)
            self._a    = 0 if astr.isspace() else int(astr)
            self._meta = zsymam_str[10] == 'M'
        except:
            if verbose:
                print(warning("Fail to convert ZSYMAM data"))
    
    def ss(self) -> str:
        return self._ss
    
    def z(self) -> int:
        return self._z
    
    def a(self) -> int:
        return self._a
    
    def meta(self) -> bool:
        return self._meta
    
    def __repr__(self) -> str:
        return self._repr


class ANGULAR_DIST_TYPE(Enum):
    ALL_ISOTROPIC = 0  # All angular distributions are isotropic
    LEGENDRE      = 1  # Data are given as Legendre expansion coefficients
    TABULATED     = 2  # Data are given as tabulated probability distributions
    MIXED         = 3  # Legendre in low energy region, tabulated in high energy region


class INTERPOLATION_SCHEME(Enum):
    CONSTANT      = 1
    LINEAR_LINEAR = 2
    LINEAR_LOG    = 3
    LOG_LINEAR    = 4
    LOG_LOG       = 5