from __future__ import annotations

"""
ENDF main class

this code is part of the RT2 project
"""

from src.endf.stream import ENDFIfstream, RecCont
from src.endf.file1 import DescData
from src.endf.file3 import CrossSection
from src.endf.file4 import AngularDist


class Reaction:
    def __init__(self, mt):
        self.mt  = mt
        self.xs  = None  # FILE3
        self.ad  = None  # FILE4
        self.ed  = None  # FILE5 
        self.ead = None  # FILE6


class ENDF:
    def __init__(self, file_name: str, verbose: bool = False, read_only_header: bool = False):
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
                    if read_only_header:
                        break
                else:
                    pass
            elif mf == 3:  # FILE3
                if mt not in self._reactions.keys():
                    self._reactions[mt] = Reaction(mt)
                self._reactions[mt].xs = CrossSection(RecCont(text), stream, mt)
            elif mf == 4:
                if mt not in self._reactions.keys():
                    self._reactions[mt] = Reaction(mt)
                self._reactions[mt].ad = AngularDist(RecCont(text), stream, mt)
            else:
                pass
        
        stream.close()

    def __getitem__(self, mt: int) -> CrossSection:
        return self._reactions[mt]
    
    def desc(self) -> DescData:
        return self._desc
    
    def keys(self) -> list:
        return self._reactions.keys()