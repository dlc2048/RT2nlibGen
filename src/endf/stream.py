from __future__ import annotations


from src.endf.desc import REACTION_TYPE
from src.endf.record import RecText, RecCont, RecList, RecTab1, RecTab2


class ENDFIfstream:
    def __init__(self, file_name: str, verbose: bool = False):
        self._verbose = verbose
        self._stream  = open(file_name, mode='r')
        self.log("Read ENDF file {}".format(file_name))

    def close(self):
        self._stream.close()

    def log(self, line: str):
        if self._verbose:
            print(line)

    def _getline(self) -> RecText | None:
        line = self._stream.readline()
        if not line:
            return None
        return RecText(line)

    def text(self) -> RecText | None:
        return self._getline()
    
    def cont(self) -> RecCont | None:
        text = self._getline()
        if not text:
            return None
        return RecCont(text)
        
    def head(self) -> RecCont | None:
        return self.cont()
    
    def end(self) -> RecCont | None:
        return self.cont()
    
    def list(self) -> RecList | None:
        text = self._getline()
        if not text:
            return None
        rec_list = RecList(text)
        npl      = rec_list.n1()
        for i in range((npl - 1) // 6 + 1):
            rec_list.append(self._getline(), npl - 6 * i)
        return rec_list
    
    def tab1(self) -> RecTab1 | None:
        text = self._getline()
        if not text:
            return None
        rec_tab1 = RecTab1(text)
        nr2      = rec_tab1.n1() * 2
        np2      = rec_tab1.n2() * 2
        for i in range((nr2 - 1) // 6 + 1):
            rec_tab1.appendInterp(self._getline(), nr2 - 6 * i)
        for i in range((np2 - 1) // 6 + 1):
            rec_tab1.appendTab(self._getline(), np2 - 6 * i)
        return rec_tab1
    
    def tab2(self, is_list: bool, tab2_nested: bool) -> RecTab2 | None:
        text = self._getline()
        if not text:
            return None
        rec_tab2 = RecTab2(text)
        nr2      = rec_tab2.n1() * 2
        nz       = rec_tab2.n2()
        for i in range((nr2 - 1) // 6 + 1):
            rec_tab2.append(self._getline(), nr2 - 6 * i)
        if tab2_nested:
            for i in range(nz):
                rec_tab2.append(self.tab2(is_list, False))
        else:
            if is_list:
                for i in range(nz):
                    rec_tab2.append(self.list())
            else:
                for i in range(nz):
                    rec_tab2.append(self.tab1())
        return rec_tab2
    

class FileInterface:
    def __init__(self, mf: int, mt: int):
        self._mt   = mt
        self._mf   = mf
        self._repr = REACTION_TYPE[mt]

    def _checkSectionEnd(self, stream: ENDFIfstream):
        if stream.end().mt() == 0:
            return
        raise ValueError("Error found in ENDF {} {} section end".format(self._mt, self._mf))
    
    def mt(self) -> int:
        return self._mt
    
    def mf(self) -> int:
        return self._mf
    
    def repr(self) -> str:
        return self._repr
    

class CrossSectionInterface(FileInterface):
    def __init__(self, head: RecCont, stream: ENDFIfstream, mt: int, mf: int):
        super().__init__(mt, mf)
        self._za   = head.c1()
        self._mass = head.c2()
        self._tab  = stream.tab1()
        self._checkSectionEnd()

    def za(self) -> float:
        return self._za
    
    def mass(self) -> float:
        return self._mass
    
    def tab(self) -> RecTab1:
        return self._tab
    