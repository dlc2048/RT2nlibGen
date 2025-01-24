from __future__ import annotations

"""
ENDF main class

this code is part of the RT2 project
"""

from src.endf.stream import ENDFIfstream, RecCont
from src.endf.file1 import DescData
from src.endf.file3 import CrossSection


class ENDF:
    def __init__(self, file_name: str, verbose: bool = False):
        stream = ENDFIfstream(file_name, verbose)

        self._reactions = {}

        while True:
            text = stream.text()
            if text is None:
                break
            mf = text.mf()
            mt = text.mt()
            if mf == 0:  # Section, material or file end
                pass
            elif mf == 1:  # FILE1
                if mt == 451:
                    self._desc = DescData(RecCont(text), stream)
                else:
                    pass
            elif mf == 3:  # FILE3
                self._reactions[mt] = CrossSection(RecCont(text), stream, mt)
            else:
                pass

    def __getitem__(self, mt: int) -> CrossSection:
        return self._reactions[mt]
    
    def desc(self) -> DescData:
        return self._desc
    
    def keys(self) -> list:
        return self._reactions.keys()