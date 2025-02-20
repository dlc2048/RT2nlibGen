from __future__ import annotations

import copy

import numpy as np
from tqdm import tqdm

from src.endf.endf import ENDF
from src.endf.stream import ENDFIfstream, FileInterface
from src.endf.record import RecText, RecCont
from src.endf.desc import REACTION_TYPE
from src.gendf.desc import GENDF_MF_TYPE, GENDF_MF_SYMBOL, GENDF_MF_TO_MULT, GENDF_MF_TO_ZA, GENDF_ZA_TO_MF

from src.constants import MAX_MULTIPLICITY, MIN_GENION_ENERGY, SCATTERING_INTERP_N
from src.algorithm import legendreToEquibin
from src.physics import muCMToLab, muCMToLabTarget
from src.prompt import info, warning

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
        self._equiprob     = None  # Equiprobable angle bins
        self._repr         = GENDF_MF_TYPE[mf]

    def control(self) -> np.ndarray:
        return self._control
    
    def matrix(self) -> np.ndarray:
        return self._matrix
    
    def equiprobable(self) -> np.ndarray:
        return self._equiprob
    
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
    
    def generateEquiprobAngles(self, nebins: int, verbose: bool):
        self._equiprob = np.empty((self._matrix.shape[0], nebins + 1))
        modifier       = (np.arange(0, self._matrix.shape[1], 1) * 2 + 1) / 2   # ENDF legendre coeff
        iterator = range(self._matrix.shape[0])
        if verbose:
            iterator = tqdm(iterator)
            iterator.set_description('{:<10}'.format(GENDF_MF_TYPE[self._mf]))
        for i in iterator:
            poly        = self._matrix[i]
            ebin_seg, _ = legendreToEquibin(poly * modifier, nebins)
            self._equiprob[i] = ebin_seg
        return
    
    def setControl(self, cont: np.ndarray):
        self._control = cont

    def setMatrix(self, mat: np.ndarray):
        self._matrix = mat

    def setEquiprobable(self, equi: np.ndarray):
        self._equiprob = equi

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
            
            mult = total / xs
            if mult > MAX_MULTIPLICITY:  # NJOY21 conversion error
                mult = 1.0
            self._multiplicity[group] = 0 if xs == 0.0 else mult
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
    __egn = None
    __egg = None

    @staticmethod
    def setNeutronGroupStructure(egn: np.ndarray):
        Reaction.__egn = np.copy(egn)
    
    @staticmethod
    def setGammaGroupStructure(egg: np.ndarray):
        Reaction.__egg = np.copy(egg)

    @staticmethod
    def ngn():
        return len(Reaction.__egn) - 1
    
    @staticmethod
    def ngg():
        return len(Reaction.__egg) - 1

    def __init__(self, za_target: int, mt: int, qvalue: float):
        self._za_origin = za_target
        self._za_res    = None
        self._repr   = REACTION_TYPE[mt]
        self._mt     = mt
        self._qvalue = qvalue  # reaction Q value
        self._inst   = None    # reaction instruction
        self._xs     = None    # XS vector
        self._depo   = None    # Multiplicity vector
        self._mult   = None    # Deposition vector
        self._comp   = {}      # Secondaries

    def setXS(self, xs_file: CommonFile):
        xs_gen = CrossSection(xs_file)
        self._xs   = xs_gen.xs()
        self._mult = np.zeros_like(self._xs)   
        self._depo = np.zeros_like(self._xs)   

    def setXSArr(self, xs: np.ndarray):
        self._xs = xs

    def addComponent(self, mf: int, trans_file: CommonFile):
        self._comp[mf] = SecondaryFactory.interpret(mf, Reaction.ngn(), Reaction.ngg(), trans_file)

    def setComponent(self, mf: int, secondary: Secondary):
        self._comp[mf] = secondary

    def residualZA(self) -> int | None:
        return self._za_res

    def _normalizeXS(self):
        for mf in self._comp:
            comps = self._comp[mf]
            for comp in comps:
                if comp is not None:
                    if self.xs() is None:  # MF13 case -> get XS from MF13
                        self._xs   = comp.normalizeMF13()
                        self._mult = np.zeros_like(self._xs)   
                        self._depo = np.zeros_like(self._xs)
                    else:  # get XS from MF3
                        comp.normalize(self.xs())
    
    def _mergeTransitionmatrix(self):
        for mf in self._comp:
            trans, spect = self._comp[mf]
            if spect is None:
                self._comp[mf] = trans
            elif trans is None:
                self._comp[mf] = spect
            else:
                trans.merge(spect)
                self._comp[mf] = trans

    def _prepareInstruction(self, verbose: bool):
        # calculate residual nucleus ZA
        za_res = self._za_origin + 1  # neutron capture
        for mf in GENDF_MF_TO_MULT.keys():
            multiplicity = GENDF_MF_TO_MULT[mf][self.mt()]
            if mf in self._comp.keys():
                za_res   -= multiplicity * GENDF_MF_TO_ZA[mf]

        # reaction instruction
        stape = []

        # hadron secondary
        for mf in GENDF_MF_TO_MULT.keys():
            multiplicity = GENDF_MF_TO_MULT[mf][self.mt()]
            if GENDF_MF_TO_ZA[mf] == za_res and multiplicity == 0:
                multiplicity = 1
            if multiplicity == 0:
                continue
            if mf in self._comp.keys():
                stape += [mf] * multiplicity

        # gamma
        if 16 in self._comp.keys():
            multiplicity = np.max(self._comp[16].multiplicity())
            multiplicity = int(np.ceil(multiplicity))
            stape += [16] * multiplicity
        
        self._inst = np.array(stape, dtype=int)

        # calculate za res again
        za_res = self._za_origin + 1  # neutron capture
        for mf in GENDF_MF_TO_MULT.keys():
            multiplicity = GENDF_MF_TO_MULT[mf][self.mt()]
            za_res      -= multiplicity * GENDF_MF_TO_ZA[mf]
        if za_res in GENDF_ZA_TO_MF:
            mf = GENDF_ZA_TO_MF[za_res]
            if mf in self._comp.keys():
                za_res = 0
        self._za_res = za_res

        inst_list = list(map(lambda x: GENDF_MF_SYMBOL[x], self._inst))
        print(info("Reaction instruction of MT={}, {} reaction".format(self._mt, REACTION_TYPE[self._mt])))
        print("Secondaries -> [{}]".format(', '.join(inst_list)))
        print("Residual ZA -> {}".format(self._za_res))

            
    def _calculateResDoseAndMultiplicity(self, verbose):
        
        egn = Reaction.__egn
        egg = Reaction.__egg

        egn_mean = (egn[1:] + egn[:-1]) * 0.5
        egg_mean = (egg[1:] + egg[:-1]) * 0.5

        # deposition
        if not self._za_res:  # target remnant not exist
            print(info('Residual target remnant not exist for MT={}, {} reaction. Energy deposition is set to 0.'.format(self._mt, REACTION_TYPE[self._mt])))
            self._depo[:] = 0.0
        else:
            if 26 in self._comp.keys():  # Tabulated depo data exist
                print(info('Energy deposition data found for MT={}, {} reaction'.format(self._mt, REACTION_TYPE[self._mt])))
                for group in range(Reaction.ngn()):
                    mask = self._comp[26].probabilityMask(group)
                    self._depo[group] = np.sum(egn_mean * mask)
            else:  # Not exist
                print(info('Energy deposition data not found for MT={}, {} reaction. Using Q-value to calculate deposition.'.format(self._mt, REACTION_TYPE[self._mt])))
                self._depo[:] = self._qvalue + egn_mean
                for group in range(Reaction.ngn()):
                    for mf_this in self._inst:
                        mask = self._comp[mf_this].probabilityMask(group)
                        if mf_this == 16:
                            self._depo[group] -= np.sum(egg_mean * mask) * self._comp[mf_this].multiplicity()[group]
                            break  # end if photon
                        else:
                            self._depo[group] -= np.sum(egn_mean * mask)

        self._depo[self._depo < 0.0] = 0.0
        self._depo[:] *= 1e-6  # eV to MeV

        # remove 26 (res-dose) if exist
        if 26 in self._comp.keys():
            del self._comp[26]

        # multiplicity
        n_hadron      = np.sum(self._inst != 16)
        multiplicity  = self._comp[16].multiplicity() if 16 in self._inst else 0.0
        multiplicity += n_hadron
        self._mult[:] = multiplicity

    def _generateElasticEquiprobAngles(self, endf: ENDF, nebins: int, len_thermal: int, verbose: bool):
        print(info('Angular distribution is generated from the scattering kinematics for MT={}, {} reaction'.format(self._mt, REACTION_TYPE[self._mt])))
        
        # generate elastic equiprob for the scattered neutron
        matrix     = self._comp[6].matrix()
        control    = self._comp[6].control()

        equiprob_n = np.empty((matrix.shape[0], nebins + 1))
        modifier   = (np.arange(0, matrix.shape[1], 1) * 2 + 1) / 2   # ENDF legendre coeff

        # thermal -> using legendre polynomial
        cpos, _, len_c = control[len_thermal]
        iterator = range(cpos)
        if verbose:
            iterator = tqdm(iterator)
            iterator.set_description('Thermal')
        for i in iterator:
            ebin_seg, _ = legendreToEquibin(matrix[i] * modifier, nebins)
            equiprob_n[i] = ebin_seg

        # fast -> kinematics
        egn   = Reaction.__egn
        awr   = endf.desc().mass()
        awr   = max(1.0, awr)

        ad_mt2 = endf[2].ad
        is_cm  = ad_mt2.isOnCMSystem()

        iterator = range(len_thermal, Reaction.ngn())
        if verbose:
            iterator = tqdm(iterator)
            iterator.set_description('Fast scattered')
        for gin in iterator:
            cpos, gmin, len_c = control[gin]
            if len_c < 0:
                continue

            # get transition prob
            pseg = matrix[cpos:cpos + len_c, 0]
            pseg = pseg / np.sum(pseg)  # normalize

            ein  = (egn[gin] + egn[gin + 1]) * 0.5
            dist = ad_mt2.getDistribution(ein)
            
            # get mu bins from transition probability (scattered neutron)
            mu_min = -1.0
            mu_max = +1.0
            if not is_cm:
                mu_min = muCMToLab(awr, mu_min)
                mu_max = muCMToLab(awr, mu_max)

            ebin_seg_list, _ = dist.getEquibin(nebins, mu_min, mu_max, pseg)
            if is_cm:
                ebin_seg_list = muCMToLab(awr, ebin_seg_list)

            equiprob_n[cpos:cpos + len_c] = ebin_seg_list
        
        # overwrite data
        self._comp[6].setMatrix(matrix)
        self._comp[6].setEquiprobable(equiprob_n)

    def _generateEvapEquiprobAngles(self, nebins: int, verbose: bool):
        print(info('Angular distribution is generated from the GENDF for MT={}, {} reaction'.format(self._mt, REACTION_TYPE[self._mt])))
        for mf in self._comp.keys():
            if mf not in (16, 26):  # only for hadron
                self._comp[mf].generateEquiprobAngles(nebins, verbose)
    
    def normalize(self):
        if self._mt == 1:  # No rule for (n,total) MT
            return
        self._normalizeXS()
        self._mergeTransitionmatrix()

    def prepareData(self, verbose: bool):
        if self._mt == 1:  # No rule for (n,total) MT
            return
        self._prepareInstruction(verbose)
        self._calculateResDoseAndMultiplicity(verbose)

    def generateEquiprobAngles(self, endf: ENDF, nebins: int, len_thermal: int, verbose: bool):
        if self._mt == 1:  # No rule for (n,total) MT
            return
        elif self._mt == 2:  # elastic
            self._generateElasticEquiprobAngles(endf, nebins, len_thermal, verbose)
        else:
            self._generateEvapEquiprobAngles(nebins, verbose)
        

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
    
    def inst(self) -> list:
        return self._inst
    
    def xs(self) -> np.ndarray | None:
        return self._xs
    
    def depo(self) -> np.ndarray | None:
        return self._depo
    
    def multiplicity(self) -> np.ndarray | None:
        return self._mult
    

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