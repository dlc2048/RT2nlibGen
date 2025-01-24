from __future__ import annotations

"""
ENDF Flags

this code is part of the RT2 project
"""


class ZSYMAM:
    def __init__(self, zsymam_str: str):
        self._repr = zsymam_str

        zstr = zsymam_str[0:3]
        astr = zsymam_str[7:10]

        self._ss   = zsymam_str[4:6]
        self._z    = 0 if zstr.isspace() else int(zstr)
        self._a    = 0 if astr.isspace() else int(astr)
        self._meta = zsymam_str[10] == 'M'
    
    def ss(self) -> str:
        return self._ss
    
    def z(self) -> int:
        return self._z
    
    def a(self) -> int:
        return self._a
    
    def meta(self) -> bool:
        return self._meta

    