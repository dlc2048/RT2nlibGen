from __future__ import annotations

"""
ENDF FILE 3: REACTION CROSS SECTIONS

this code is part of the RT2 project
"""

from src.endf.stream import CrossSectionInterface, ENDFIfstream
from src.endf.record import RecCont


class CrossSection(CrossSectionInterface):
    def __init__(self, head: RecCont, stream: ENDFIfstream, mt: int):
        super().__init__(head, stream, mt, 3)

    def massQvalue(self) -> float:
        return self._tab.c1()
    
    def reactionQvalue(self) -> float:
        return self._tab.c2()
    
    def isBreakup(self) -> bool:
        return bool(self._tab.l2())