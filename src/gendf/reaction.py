from __future__ import annotations

import copy

import numpy as np

from src.endf.stream import ENDFIfstream, FileInterface
from src.endf.record import RecText, RecCont
from src.endf.desc import REACTION_TYPE
from src.gendf.desc import GENDF_MF_TYPE, GENDF_MF_TO_MULT

class CommonFile(FileInterface):
    def __init__(self, head: RecText, stream: ENDFIfstream):
        super().__init__(head.mf(), head.mt())
        cont         = RecCont(head)
        self._za     = int(cont.c1())
        self._awr    = cont.c2()
        self._nl     = int(cont.l1())
        self._nz     = int(cont.l2())
        self._lrflag = bool(cont.n1())
        self._ngn    = cont.n2()
        self._list   = []
        while True:
            self._list += [stream.list()]
            if self._list[-1].n2() == self._ngn:
                break
        
        self._checkSectionEnd(stream)

    def legendreOrder(self) -> int:
        return self._nl
    
    def ngroup(self) -> int:
        return self._ngn

    def items(self) -> list:
        return self._list
    
    def setItems(self, items: list):
        self._list = items


class CrossSection:
    def __init__(self, xs_file: CommonFile):
        self._repr   = GENDF_MF_TYPE[3]
        self._ngroup = xs_file.ngroup()
        self._xs     = np.zeros(self._ngroup, dtype=float)
        for seg in xs_file.items():
            group = seg.n2() - 1
            sig   = seg.list()[1]
            self._xs[group] = sig

    def xs(self) -> np.ndarray:
        return self._xs
    
    def ngroup(self) -> float:
        return self._ngroup
    
    def __repr__(self) -> str:
        return self._repr
    

class Secondary:
    def __init__(self, mf: int, ngn: int, ngg: int):
        self._mf           = mf
        self._ngn          = ngn
        self._ngg          = ngg
        self._multiplicity = np.zeros(ngn)
        self._control      = -np.ones((ngn, 3), dtype=int)     # memory control  [matrix position, lowest transition group, matrix length]
        self._matrix       = None  # Legendre polynomials
        self._repr         = GENDF_MF_TYPE[mf]

    def control(self) -> np.ndarray:
        return self._control
    
    def matrix(self) -> np.ndarray:
        return self._matrix
    
    def multiplicity(self) -> np.ndarray:
        return self._multiplicity
    
    def probabilityMask(self, group: int) -> np.ndarray:
        mask_dim = self._ngn if self._mf != 16 else self._ngg
        mask     = np.zeros(mask_dim)
        pos, lgroup, length = self._control[group]
        if length >= 0:
            mask[lgroup:lgroup + length] = self._matrix[pos:pos + length, 0]
            mask /= np.sum(mask)  # normalize
        return mask
    
    def setControl(self, cont: np.ndarray):
        self._control = cont

    def setMatrix(self, mat: np.ndarray):
        self._matrix = mat

    def setMultiplicity(self, multiplicity: np.ndarray):
        self._multiplicity = multiplicity


class Transition(Secondary):
    def __init__(self, mf: int, ngn: int, ngg: int, trans_file: CommonFile):
        super().__init__(mf, ngn, ngg)
        nl     = trans_file.legendreOrder()
        ntrans = 0

        for seg in trans_file.items():
            ntrans += int(seg.l1()) - 1

        self._matrix = np.empty((ntrans, nl), dtype=float)  # Legendre polynomials
        cpos = 0  # current matrix position
        for seg in trans_file.items():
            group   = seg.n2() - 1
            ig2lo   = int(seg.l2()) - 1
            mat_seg = seg.list()
            mat_seg = mat_seg.reshape((-1, nl))[1:]  # 1D to 2D
            mat_len = mat_seg.shape[0]
            self._control[group] = (cpos, ig2lo, mat_len)
            self._matrix[cpos:cpos + mat_len] = mat_seg
            cpos += mat_len

    def normalize(self, xs_arr: np.ndarray):
        # normalize transition sampling table
        mask = np.zeros(self._matrix.shape[0], dtype=bool)
        for group, control in enumerate(self._control):
            cpos, _, mat_len = control
            xs = xs_arr[group]
            if cpos < 0:
                continue
            if np.any(mask[cpos:cpos + mat_len]):
                raise ValueError('Transition sampling table is corrupted in {}'.format(self.__repr__()))
            mask[cpos:cpos + mat_len] = True
            seg      = self._matrix[cpos:cpos + mat_len]
            leading  = seg[:,0:1]
            total    = np.sum(leading)

            if total == 0.0:  # Invalid group 
                self._control[group] = (-1, -1, -1)
                continue
                
            self._multiplicity[group] = 0 if xs == 0.0 else total / xs
            seg /= total

    def normalizeMF13(self) -> np.ndarray:
        mask   = np.zeros(self._matrix.shape[0], dtype=bool)
        ngroup = self._control.shape[0]
        xs     = np.zeros(ngroup, dtype=float)
        for group, control in enumerate(self._control):
            cpos, _, mat_len = control
            if cpos < 0:
                continue
            if np.any(mask[cpos:cpos + mat_len]):
                raise ValueError('Transition sampling table is corrupted in {}'.format(self.__repr__()))
            mask[cpos:cpos + mat_len] = True
            seg       = self._matrix[cpos:cpos + mat_len]
            leading   = seg[:,0:1]
            total     = np.sum(leading)
            xs[group] = total

            if total == 0.0:  # Invalid group 
                self._control[group] = (-1, -1, -1)
                continue

            self._multiplicity[group] = 0 if total == 0.0 else 1.0
            seg /= total
        return xs

    def merge(self, spect: Spectrum):
        mat1  = spect.matrix()
        cont1 = spect.control()
        mul1  = spect.multiplicity()
        mat2  = self.matrix()
        cont2 = self.control()
        mul2  = self.multiplicity()
        dim_max = max(mat1.shape[1], mat2.shape[1])
        offset  = mat1.shape[0]

        mat_merged  = np.zeros((mat1.shape[0] + mat2.shape[0], dim_max))
        mat_merged[:offset,:mat1.shape[1]] = mat1
        mat_merged[offset:,:mat2.shape[1]] = mat2

        mask1 = cont1[:,0] >= 0
        mask2 = cont2[:,0] >= 0

        # control card
        cont_merged = np.ones_like(cont1)
        cont_merged[mask1]    = cont1[mask1]
        cont_merged[mask2]    = cont2[mask2]
        cont_merged[mask2,0] += offset
        # check integrity
        if np.any(mask1 * mask2):
            raise ValueError('Merging failed in {}. Conflicted group'.format(self.__repr__()))
        
        # multiplicity
        mul_merged = np.zeros_like(mul1)
        mul_merged[mask1] = mul1[mask1]
        mul_merged[mask2] = mul2[mask2]

        self._control      = cont_merged
        self._matrix       = mat_merged
        self._multiplicity = mul_merged

        return

    def __repr__(self) -> str:
        return '{} transition'.format(self._repr)
    

class Spectrum(Secondary):
    def __init__(self, mf: int, ngn: int, ngg: int, spect_file: CommonFile):
        super().__init__(mf, ngn, ngg)
        nl      = spect_file.legendreOrder()
        control = None
        for seg in spect_file.items():
            ig2lo = int(seg.l2())
            ig    = seg.n2()
            if not ig:  # spectrum
                spect_len    = seg.l1()
                self._matrix = np.reshape(seg.list(), (spect_len, nl))  # set spectrum matrix
                control      = (0, ig2lo - 1, spect_len)
            elif not ig2lo:  # multiplicity
                group     = seg.n2() - 1
                multi_seg = seg.list()[1]
                self._multiplicity[group] = multi_seg

        # fill control
        for group in range(ngn):
            if self._multiplicity[group] > 0:
                self._control[group] = control
        return
    
    def normalize(self, xs: np.ndarray):
        # normalize multiplicity
        self._multiplicity = np.divide(self._multiplicity, xs, out=np.zeros_like(xs),where=xs!=0)
        # normalize spectrum sampling table
        mask = np.zeros(self._matrix.shape[0], dtype=bool)
        for cpos, _, mat_len in self._control:
            if cpos < 0:
                continue
            if np.all(mask[cpos:cpos + mat_len]):
                continue
            if np.any(mask[cpos:cpos + mat_len]):
                raise ValueError('Spectrum sampling table is corrupted in {}'.format(self.__repr__()))
            mask[cpos:cpos + mat_len] = True
            seg      = self._matrix[cpos:cpos + mat_len]
            leading  = seg[:,0:1]
            seg     /= np.sum(leading)
        return

    def __repr__(self) -> str:
        return  '{} spectrum'.format(self._repr)


class SecondaryFactory:
    @staticmethod
    def interpret(mf: int, ngn: int, ngg: int, common_file: CommonFile):
        items = common_file.items()
        items_trans = []
        items_spect = []
        while len(items):
            seg = items.pop()
            ig2lo = int(seg.l2())
            ig    = seg.n2()
            if not ig2lo or not ig:
                items_spect = [seg] + items_spect
            else: 
                items_trans = [seg] + items_trans
        
        # rebuild
        trans_file = copy.deepcopy(common_file)
        spect_file = copy.deepcopy(common_file)
        trans_file.setItems(items_trans)
        spect_file.setItems(items_spect)

        transition = Transition(mf, ngn, ngg, trans_file) if len(items_trans) else None
        spectrum   = Spectrum(mf, ngn, ngg, spect_file)   if len(items_spect) else None
        return transition, spectrum


class Reaction:
    def __init__(self, mt: int, ngn: int, ngg: int, qvalue: float):
        self._repr   = REACTION_TYPE[mt]
        self._mt     = mt
        self._ngn    = ngn
        self._ngg    = ngg
        self._qvalue = qvalue
        self._xs     = None
        self._mult   = None
        self._comp   = {}

    def setXS(self, xs_file: CommonFile):
        xs_gen = CrossSection(xs_file)
        self._xs   = xs_gen.xs()               # XS vector
        self._mult = np.zeros(self._xs.shape)  # Multiplicity vector

    def addComponent(self, mf: int, trans_file: CommonFile):
        self._comp[mf] = SecondaryFactory.interpret(mf, self._ngn, self._ngg, trans_file)

    def setComponent(self, mf: int, secondary: Secondary):
        self._comp[mf] = secondary

    def normalize(self):
        for mf in self._comp:
            comps = self._comp[mf]
            for comp in comps:
                if comp is not None:
                    if self.xs() is None:  # MF13 case -> get XS from MF13
                        self._xs = comp.normalizeMF13()
                    else:  # get XS from MF3
                        comp.normalize(self.xs())
    
    def merge(self):
        for mf in self._comp:
            trans, spect = self._comp[mf]
            if spect is None:
                self._comp[mf] = trans
            elif trans is None:
                self._comp[mf] = spect
            else:
                trans.merge(spect)
                self._comp[mf] = trans

    def __repr__(self) -> str:
        return self._repr
    
    def __getitem__(self, index):
        return self._comp[index]
    
    def keys(self):
        return self._comp.keys()
    
    def mt(self) -> int:
        return self._mt
    
    def qvalue(self) -> float:
        return self._qvalue
    
    def xs(self) -> np.ndarray | None:
        return self._xs
    
    def ngn(self) -> int:
        return self._ngn
    
    def ngg(self) -> int:
        return self._ngg
    

def mergeComponents(mf: int, reactions: list) -> Secondary:
    ngn = reactions[0].ngn()
    ngg = reactions[0].ngg()

    cont     = -np.ones((ngn, 3), dtype=int)
    mat_list = []
    mult     = np.zeros(ngn, dtype=float)

    # check legendre & total XS
    max_lorder = -1
    xs_total   = np.zeros(ngn, dtype=float)
    for reaction in reactions:
        xs_total += reaction.xs()
        if mf not in reaction.keys():
            continue
        max_lorder = max(max_lorder, reaction[mf].matrix().shape[1])

    mat_len_total = 0
    for group in range(ngn):
        # check gmin (ig2lo) and gmax
        gmin = 999999999
        gmax = -1
        for reaction in reactions:
            if mf not in reaction.keys():
                continue
            cpos, ig2lo, mat_len = reaction[mf].control()[group]
            if mat_len < 0:
                continue
            gmin = min(gmin, ig2lo)
            gmax = max(gmax, ig2lo + mat_len)
        if gmax - gmin <= 0:
            continue
        # control
        mat_seg        = np.zeros((gmax - gmin, max_lorder), dtype=float)
        cont[group]    = [mat_len_total, gmin, gmax - gmin]
        mat_len_total += gmax - gmin
        # matrix
        for reaction in reactions:
            if mf not in reaction.keys():
                continue
            cpos, ig2lo, mat_len = reaction[mf].control()[group]
            if mat_len < 0:
                continue
            xs    = reaction.xs()[group]
            cmult = 0
            if mf == 16:  # gamma
                cmult = reaction[mf].multiplicity()[group]
                if cmult == 0.0:
                    cmult = 1.0
            elif mf == 26:  # depo
                cmult = 1.0
            else:
                cmult = GENDF_MF_TO_MULT[mf][reaction.mt()]
            mat_sub      = reaction[mf].matrix()
            mult[group] += cmult * xs / xs_total[group]
            for i in range(mat_len):
                poly = mat_sub[cpos + i]
                gpos = ig2lo + i - gmin
                mat_seg[gpos,:len(poly)] += poly * xs * cmult
        mat_list += [mat_seg]

    mat    = np.zeros((mat_len_total, max_lorder), dtype=float)
    offset = 0
    for mat_seg in mat_list:
        mat[offset:offset + mat_seg.shape[0]] = mat_seg[:]
        offset += mat_seg.shape[0]
    
    secondary = Secondary(mf, ngn, ngg)
    secondary.setMultiplicity(mult)
    secondary.setControl(cont)
    secondary.setMatrix(mat)
    return secondary