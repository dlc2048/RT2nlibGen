from __future__ import annotations

"""
ENDF record elements

this code is part of the RT2 project
"""

import numpy as np


def stringToEndfFloat(string: str) -> float:
    """
    convert endf-float format to scientific notation
    """
    if "Inf" in string:
        string = "0.0"
    if not string.strip():
        string = "0.0"
    if "." in string:
        string = string[0] + string[1:].replace("+", "e+", 1)
        string = string[0] + string[1:].replace("-", "e-", 1)
    return float(string)


def stringToEndfInt(string: str) -> int:
    if not string.strip():
        return 0
    else:
        return int(string)


class RecordInterface:
    def __init__(self, line: str):
        self._mat = int(line[66:70])
        self._mf  = int(line[70:72])
        self._mt  = int(line[72:75])
        self._ns  = -1 if line[75:80].isspace() else int(line[75:80])

    def mat(self) -> int:
        return self._mat
    
    def mf(self) -> int:
        return self._mf
    
    def mt(self) -> int:
        return self._mt
    
    def ns(self) -> int:
        return self._ns
    

class RecText(RecordInterface):
    def __init__(self, line: str):
        super().__init__(line)
        self._hl = line[0:66]

    def hl(self) -> str:
        return self._hl
    

class RecCont(RecordInterface):
    def __init__(self, rec_text: RecText):
        self._mat = rec_text.mat()
        self._mf  = rec_text.mf()
        self._mt  = rec_text.mt()
        self._ns  = rec_text.ns()

        contents  = rec_text.hl()
        self._c1  = stringToEndfFloat(contents[00:11])
        self._c2  = stringToEndfFloat(contents[11:22])
        self._l1  = stringToEndfInt(contents[22:33])
        self._l2  = stringToEndfInt(contents[33:44])
        self._n1  = stringToEndfInt(contents[44:55])
        self._n2  = stringToEndfInt(contents[55:66])

    def c1(self) -> float:
        return self._c1
    
    def c2(self) -> float:
        return self._c2
    
    def l1(self) -> int:
        return self._l1
    
    def l2(self) -> int:
        return self._l2
    
    def n1(self) -> int:
        return self._n1
    
    def n2(self) -> int:
        return self._n2
    

class RecList(RecCont):
    def __init__(self, rec_text: RecText):
        super().__init__(rec_text)
        self._list = np.empty(self.n1())
        self._pos  = 0
    
    def append(self, rec_text: RecText, max_inst: int = -1):
        if max_inst < 0 or max_inst > 6:
            max_inst = 6
        contents = rec_text.hl()
        for i in range(max_inst):
            self._list[self._pos] = stringToEndfFloat(contents[i * 11:(i + 1) * 11])
            self._pos += 1

    def list(self) -> np.ndarray:
        return self._list


class RecTab1(RecCont):
    def __init__(self, rec_text: RecText):
        super().__init__(rec_text)
        self._interp = np.empty((self.n1(), 2), dtype=int)
        self._tab    = np.empty((self.n2(), 2), dtype=float)
        self._ipos   = 0
        self._tpos   = 0

    def appendInterp(self, rec_text: RecText, max_inst: int = -1):
        if max_inst < 0 or max_inst > 6:
            max_inst = 6
        contents = rec_text.hl()
        for i in range(0, max_inst, 2):
            for j in range(2):
                self._interp[self._ipos, j] = stringToEndfInt(contents[(i + j) * 11:(i + j + 1) * 11])
            self._ipos += 1

    def appendTab(self, rec_text: RecText, max_inst: int = -1):
        if max_inst < 0 or max_inst > 6:
            max_inst = 6
        contents = rec_text.hl()
        for i in range(0, max_inst, 2):
            for j in range(2):
                self._tab[self._tpos, j] = stringToEndfFloat(contents[(i + j) * 11:(i + j + 1) * 11])
            self._tpos += 1

    def interp(self) -> np.ndarray:
        return self._interp
    
    def tab(self) -> np.ndarray:
        return self._tab
    

class RecTab2(RecCont):
    def __init__(self, rec_text: RecText):
        super().__init__(rec_text)
        self._interp = np.empty((self.n1(), 2), dtype=int)
        self._ipos   = 0
        self._tab    = []

    def appendInterp(self, rec_text: RecText, max_inst: int = -1):
        if max_inst < 0 or max_inst > 6:
            max_inst = 6
        contents = rec_text.hl()
        for i in range(0, max_inst, 2):
            for j in range(2):
                self._interp[self._ipos, j] = stringToEndfInt(contents[(i + j) * 11:(i + j + 1) * 11])
            self._ipos += 1

    def append(self, tab: RecTab2 | RecTab1 | RecList):
        self._tab += [tab]

    def interp(self) -> np.ndarray:
        return self._interp
    
    def tab(self) -> list:
        return self._tab

    def good(self) -> bool:
        return len(self.list()) ^ len(self.tab1())