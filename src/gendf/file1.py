from __future__ import annotations

import numpy as np

from src.endf.stream import ENDFIfstream, FileInterface
from src.endf.record import RecText, RecCont


class DescData(FileInterface):
    def __init__(self, head: RecText, stream: ENDFIfstream):
        super().__init__(1, 451)
        cont      = RecCont(head)
        self._za  = int(cont.c1())
        self._awr = cont.c2()
        self._nz  = cont.l2()
        self._ntw = cont.n2()
        f1_list   = stream.list()
        ebins     = f1_list.list()
        offset    = self._nz + self._ntw
        self._ngn = f1_list.l1()
        self._ngg = f1_list.l2()
        self._egn = ebins[offset:offset + self._ngn + 1]
        self._egg = ebins[offset + self._ngn + 1:]
        self._checkSectionEnd(stream)

    def ngn(self) -> int:
        return self._ngn
    
    def ngg(self) -> int:
        return self._ngg
    
    def egn(self) -> np.ndarray:
        return self._egn
    
    def egg(self) -> np.ndarray:
        return self._egg
    
    def za(self) -> int:
        return self._za
    