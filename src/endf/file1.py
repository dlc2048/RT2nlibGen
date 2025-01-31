from __future__ import annotations

"""
ENDF FILE 1: GENERAL INFORMATION

this code is part of the RT2 project
"""


from src.endf.stream import ENDFIfstream, FileInterface
from src.endf.record import RecCont
from src.endf.param import ZSYMAM
from src.endf.desc import REACTION_TYPE


class DescData(FileInterface):
    def __init__(self, head: RecCont, stream: ENDFIfstream, verbose: bool = False):
        super().__init__(1, 451)
        self._mat  = head.mat()
        self._za   = int(head.c1())
        self._awr  = head.c2()
        self._lrp  = head.l1()
        self._lfi  = bool(head.l2())
        self._nlib = head.n1()
        self._nmod = head.n2()
        # start the streaming
        cont = stream.cont()
        self._elis = cont.c1()
        self._sta  = bool(cont.c2())
        self._lis  = cont.l1()
        self._liso = cont.l2()
        self._nfor = cont.n2()
        # Third line
        cont = stream.cont()
        self._awi  = cont.c1()
        self._emax = cont.c2()
        self._lrel = cont.l1()
        self._nsub = cont.n1()
        self._nver = cont.n2()
        # Fourth line
        cont = stream.cont()
        self._temp = cont.c1()
        self._ldrv = cont.l1()
        nwd = cont.n1()
        nxc = cont.n2()
        # Text
        text = stream.text()
        self._zsymam = ZSYMAM(text.hl()[:11], verbose=verbose)
        self._alab   = text.hl()[11:22]
        self._edate  = text.hl()[22:32]
        self._auth   = text.hl()[33:66]
        text = stream.text()
        self._ref    = text.hl()[1:22]
        self._ddate  = text.hl()[22:32]
        self._rdate  = text.hl()[33:43]
        self._endate = text.hl()[55:63]
        # HSUB entry
        ss = ""
        for i in range(2, nwd):
            text = stream.text()
            ss + text.hl() + "\n"
        self._hsub = ss
        # Reaction entry
        self._rlist = []
        for i in range(nxc):
            cont = stream.cont()
            self._rlist += [[cont.l1(), cont.l2(), cont.n1(), cont.n2()]]
        self._checkSectionEnd(stream)

    def za(self) -> int:
        return self._za
    
    def mat(self) -> int:
        return self._mat
    
    def mass(self) -> int:
        return self._awr
    
    def zsymam(self) -> ZSYMAM:
        return self._zsymam
    
    def fissile(self) -> bool:
        return self._lfi
    
    def unstable(self) -> bool:
        return self._sta
    
    def isomericNumber(self) -> int:
        return self._liso
    
    def emax(self) -> float:
        return self._emax
    
    def temperature(self) -> float:
        return self._temp