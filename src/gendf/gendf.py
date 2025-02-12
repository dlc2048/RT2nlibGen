from __future__ import annotations

"""
convert GENDF file that obtained from NJOY21 to RT2 groupwised neutron library

this code is part of the RT2 project
"""

import numpy as np
from tqdm import tqdm

from src.endf.endf import ENDF
from src.algorithm import AliasTable, legendreToEquibin
from src.fortran import Fortran
from src.prompt import error, info, warning
from src.endf.stream import ENDFIfstream
from src.endf.desc import REACTION_TYPE, MTHierarchy
from src.gendf.desc import GENDF_MF_TYPE, GENDF_MF_TO_MULT, GENDF_MF_TO_ZA
from src.gendf.file1 import DescData
from src.gendf.reaction import CommonFile, Reaction, mergeComponents

from src.physics import energyToMuCM, muCMToLab, muLabToCM


class GENDF:
    __egn = None

    def __init__(self, file_name: str, nebins: int, temperature: float, endf: ENDF, sab: str = "", mfactor: float = 1.0, verbose: bool = False):
        self._desc_origin = endf.desc()
        self._temperature = temperature
        self._sab         = sab + "\0"
        self._reaction    = {}
        # reaction
        self._mt_list   = []    # possible MT reaction lists
        self._ralias    = None  # Reaction alias index
        self._rprob     = None  # Reaction alias probability
        self._stape_len = []    # Length of sampling instruction
        self._stape_off = []    # Offset of sampling instruction
        self._stape     = []    # Sampling instruction
        # group
        self._multiplicity = None  # Global particle multiplicity
        self._gcontrol     = None  # Group control card
        self._edepo        = None  # Global energy deposition
        self._galias       = None  # Group alias index
        self._gprob        = None  # Group alias probability
        self._eabin        = None  # Equiprobable angle bin

        # thermal
        self._len_thermal  = None
        
        # raw data
        self._readStream(file_name, endf, verbose)
        self._mergeSpectrum()
        self._splitRecoil(verbose)
        
        # merge termal neutron scattering matrix to elastic scattering {MT=221 -> MT=2}
        if 221 in self._reaction.keys():
            self._mergeThermal(verbose, mfactor)
        # remove production cross section
        for mt in range(202, 208):
            if mt in self._reaction.keys():
                del self._reaction[mt]
        # remove anything
        if 5 in self._reaction.keys():
            del self._reaction[5]
        # merge inelastic reactions
        for mt in (4, 16, 103, 104, 105, 106, 107):
            self._mergeReaction(mt, verbose)
        # prepare reaction sampling table
        self._prepareReactionTable(verbose)
        self._prepareReactionInstruction(verbose)
        # prepare reaction offset data and empty alias & equiprob angle bin, multiplicity and energy deposition
        self._prepareUnifiedNDL(verbose)
        self._prepareUnifiedMultiplicityAndDeposition(verbose)
        # generate group alias table & equiprobable angle bins
        self._generateEquiprobAngles(endf, nebins, verbose)
        self._generateGroupAlias(verbose)

    @staticmethod
    def setNeutronGroup(egn: np.ndarray):
        GENDF.__egn = egn

    @staticmethod
    def _streamVerbose(mt: int, mf: int, verbose: bool):
        if verbose:
            print(info('Data block found for MT={}, MT={}, {} data for {}'.format(mt, mf, GENDF_MF_TYPE[mf], REACTION_TYPE[mt])))

    def _readStream(self, file_name: str, endf: ENDF, verbose: bool):
        stream = ENDFIfstream(file_name, verbose)
        while True:  # Read stream
            text = stream.text()
            if text is None:
                break
            
            mf = text.mf()
            mt = text.mt()

            if mf == 0:
                continue
            elif mf == 1:
                self._desc = DescData(text, stream)
            else:  # common file
                if mt not in self._reaction.keys():  # push new
                    q_value = 0.0
                    if mt in endf.keys():
                        q_value = endf[mt].xs.reactionQvalue()
                    
                    self._reaction[mt] = Reaction(mt, self._desc.ngn(), self._desc.ngg(), q_value)
                segs = CommonFile(text, stream)
                if mf == 3:  # XS
                    self._reaction[mt].setXS(segs)
                else:  # components
                    self._reaction[mt].addComponent(mf, segs)
            self._streamVerbose(mt, mf, verbose)
        stream.close()

    def _mergeSpectrum(self):
        for mt in self._reaction:
            self._reaction[mt].normalize()
            self._reaction[mt].merge()

    def _splitRecoil(self, verbose: bool):
        za = self._desc_origin.za()
        if za not in GENDF_MF_TO_ZA.values():
            return
        za
        pass

    def _mergeThermal(self, verbose: bool, mfactor: float):
        if verbose:
            print(info('Merge MT=221 (thermal) data to MT=2 (elastic)'))

        # append XS of MT=221 to MT=2
        thermal_xs = self._reaction[221].xs()
        elastic_xs = self._reaction[2].xs()

        len_thermal = np.argmax(thermal_xs == 0.0)

        self._len_thermal = len_thermal

        if verbose:
            print(info('Molecular number is estimated to be {}').format(mfactor))

        xs_modifier = np.copy(elastic_xs)

        elastic_xs[:len_thermal] = thermal_xs[:len_thermal] * mfactor
        self._reaction[2].setXSArr(elastic_xs)

        xs_modifier = self._reaction[2].xs() - xs_modifier

        xs_total = self._reaction[1].xs() + xs_modifier
        self._reaction[1].setXSArr(xs_total)

        thermal = self._reaction[221][6]
        elastic = self._reaction[2][6]

        cont1 = thermal.control()
        mat1  = thermal.matrix()
        mul1  = thermal.multiplicity()
        cont2 = elastic.control()
        mat2  = elastic.matrix()
        mul2  = elastic.multiplicity()

        mask1 = cont1[:,0] >= 0
        mask2 = ~mask1

        offset1 = mat1.shape[0]
        offset2 = cont2[mask2][0,0]

        # control
        cont_merge = -np.ones_like(cont1)
        cont_merge[mask1]    = cont1[mask1]
        cont_merge[mask2]    = cont2[mask2]
        cont_merge[mask2,0] += offset1 - offset2
        # check integrity
        assert np.all(mask1 + mask2) and not np.any(mask1 * mask2), error('Thermal scattering matrix summation failed in {}.'.format(self.__repr__()))
        mat_merge = np.append(mat1[:offset1], mat2[offset2:], axis=0)
        mul_merge = np.zeros_like(mul1)
        mul_merge[mask1] = mul1[mask1]
        mul_merge[mask2] = mul2[mask2]

        elastic.setControl(cont_merge)
        elastic.setMatrix(mat_merge)
        elastic.setMultiplicity(mul_merge)

        del self._reaction[221]

        return
    
    def _mergeReaction(self, mt_parent: int, verbose: bool):
        mt_list     = self._reaction.keys()
        target_list = []
        for mt in MTHierarchy[mt_parent].child():
            if mt in mt_list:
                target_list += [mt]
        if not len(target_list):
            return
        
        if verbose:
            print(info("Merge sub reactions to MT={}, {} reaction".format(mt_parent, REACTION_TYPE[mt_parent])))
            print("List of valid sub reactions")
            reaction_repr_list = ""
            for mt in target_list:
                reaction_repr_list += '{}, '.format(REACTION_TYPE[mt])
            for i in range(int(np.ceil(len(reaction_repr_list) / 80))):
                print(info('{}'.format(reaction_repr_list[i * 80: (i + 1) * 80])))

        for mf in GENDF_MF_TO_MULT:
            if GENDF_MF_TO_MULT[mf][mt] > 0:
                break

        # hadron
        parent_has_table = mf in self._reaction[mt_parent].keys()
        child_has_table  = np.any([mf in self._reaction[i].keys() for i in target_list])
        sub_reactions    = [self._reaction[i] for i in target_list]
        if not parent_has_table and child_has_table:
            secondary = mergeComponents(mf, sub_reactions)
            self._reaction[mt_parent].setComponent(mf, secondary)

        # gamma
        parent_has_table = 16 in self._reaction[mt_parent].keys()
        child_has_table  = np.any([16 in self._reaction[i].keys() for i in target_list])
        if not parent_has_table and child_has_table:
            secondary = mergeComponents(16, sub_reactions)
            self._reaction[mt_parent].setComponent(16, secondary)

        # res-dose
        parent_has_table = 26 in self._reaction[mt_parent].keys()
        child_has_table  = np.any([26 in self._reaction[i].keys() for i in target_list])
        if not parent_has_table and child_has_table:
            secondary = mergeComponents(26, sub_reactions)
            self._reaction[mt_parent].setComponent(26, secondary)

        for mt in target_list:
            del self._reaction[mt]
    
    def _prepareReactionTable(self, verbose: bool):
        # prepare target MT list
        if verbose:
            print(info("Search possible reaction branch ..."))
        mt_list = self.keys()
        mt_val  = set(mt_list)
        for mt in mt_list:
            node = MTHierarchy[mt]
            if node.parent() in mt_val:
                mt_val.remove(node.parent())

        self._mt_list = list(mt_val)
        self._mt_list.sort()

        if verbose:
            print("Possible reaction list")
            reaction_repr_list = ""
            for mt in self._mt_list:
                reaction_repr_list += '{}, '.format(REACTION_TYPE[mt])
            for i in range(int(np.ceil(len(reaction_repr_list) / 80))):
                print('{}'.format(reaction_repr_list[i * 80: (i + 1) * 80]))

        # prepare alias table
        if verbose:
            print(info("Prepare MT sampling alias table ..."))
        alias_shape  = (self._desc.ngn(), len(self._mt_list))
        self._ralias = -np.ones(alias_shape, dtype=int)
        self._rprob  = -np.ones(alias_shape, dtype=float)
        # raw probability
        prob_raw     = np.empty(alias_shape, dtype=float)
        for i, mt in enumerate(self._mt_list):
            xs = self[mt].xs()
            prob_raw[:,i] = xs
        # alias
        iterator = range(self._desc.ngn())
        if verbose:
            iterator = tqdm(iterator)
        for group in iterator:
            table = AliasTable(np.arange(len(self._mt_list)), prob_raw[group])
            self._ralias[group] = table.alias()
            self._rprob[group]  = table.prob()

    def _prepareReactionInstruction(self, verbose: bool):
        if verbose:
            print(info("Prepare reaction instructions ..."))
        for mt in self._mt_list:
            stape = []

            # hadron secondary
            for mf in GENDF_MF_TO_MULT.keys():
                multiplicity = GENDF_MF_TO_MULT[mf][mt]
                if multiplicity == 0:
                    continue
                hadron_found = False
                current_mt = mt
                while current_mt is not None:
                    if current_mt not in self._reaction.keys():
                        current_mt = MTHierarchy[current_mt].parent()
                        continue
                    if mf in self._reaction[current_mt].keys():
                        hadron_found = True
                        break
                    else:
                        current_mt = MTHierarchy[current_mt].parent()
                if hadron_found:
                    stape += [mf] * multiplicity
                else:
                    if verbose:
                        print(warning("secondary data missing for MT={}, MF={}, {} production in {}".format(mt, mf, GENDF_MF_TYPE[mf], REACTION_TYPE[mt])))

            # gamma
            multiplicity = 0.0
            gamma_found  = False
            current_mt   = mt
            while current_mt is not None:
                if current_mt not in self._reaction.keys():
                    current_mt = MTHierarchy[current_mt].parent()
                    continue
                if 16 in self._reaction[current_mt].keys():
                    gamma_found  = True
                    reaction     = self._reaction[current_mt]
                    multiplicity = np.max(reaction[16].multiplicity())
                    break
                else:
                    current_mt = MTHierarchy[current_mt].parent()
            multiplicity = int(np.ceil(multiplicity))
            stape += [16] * multiplicity if gamma_found else []  # gamma
            
            # sampling instruction
            self._stape_len += [len(stape)]
            self._stape_off += [len(self._stape)]
            self._stape     += stape

    def _prepareUnifiedNDL(self, verbose: bool):
        # reaction table -> np array
        self._mt_list   = np.array(self._mt_list, dtype=int)
        self._stape_len = np.array(self._stape_len, dtype=int)
        self._stape_off = np.array(self._stape_off, dtype=int)
        self._stape     = np.array(self._stape, dtype=int)
        # key is 1000 * mt + mf
        ntotal     = 0
        max_legend = 0  # Maximum order of legendre polynomial
        offset     = {}  # Gamma
        ngn        = self._desc.ngn()  # Number of neutron group
        # Hadron key
        for mt in self.keys():  
            reaction = self[mt]
            for mf in reaction.keys():
                mat  = reaction[mf].matrix()
                okey = 1000 * mt + mf
                if mf != 16:
                    offset[okey] = ntotal
                    ntotal      += mat.shape[0]
                    max_legend   = max(max_legend, mat.shape[1])

        self._eabin    = np.zeros((ntotal, max_legend), dtype=float)

        # Gamma key
        for mt in self.keys():  
            reaction = self[mt]
            for mf in reaction.keys():
                mat  = reaction[mf].matrix()
                okey = 1000 * mt + mf
                if mf == 16:
                    offset[okey] = ntotal
                    ntotal      += mat.shape[0]
        offset_keys = list(offset.keys())

        
        self._gcontrol = np.empty((len(offset) * ngn, 3), dtype=int)
        self._galias   = -np.ones(ntotal, dtype=int)
        self._gprob    = np.empty(ntotal, dtype=float)

        if verbose:
            print(info('Prepare unified data library ...'))
            print('Number of neutron group     : {}'.format(ngn))
            print('Length of alias table       : {}'.format(self._galias.shape[0]))
            print('Length of equiprobable table: {}'.format(self._eabin.shape[0]))
            print('Length of control table     : {}'.format(len(offset) * ngn))

        # fill control card and eabin
        if verbose:
            print(info('Prepare unified control card and polynomials ...'))
        iterator = enumerate(offset)
        if verbose:
            iterator = tqdm(iterator)
        for i, okey in iterator:
            mt = okey // 1000
            mf = okey %  1000
            reaction = self[mt]
            toff     = offset[okey]
            # control
            self._gcontrol[i * 260: (i + 1) * 260]    = reaction[mf].control()
            self._gcontrol[i * 260: (i + 1) * 260,0] += toff
            mat = reaction[mf].matrix()
            self._gprob[toff:toff + mat.shape[0]] = mat[:,0]
            if mf != 16:  # Hadron
                self._eabin[toff:toff + mat.shape[0],:mat.shape[1]] = mat

        # Link sampling instruction tape to control card
        if verbose:
            print(info('Link sampling instruction tape to control card ...'))
        iterator = range(len(self._mt_list))
        if verbose:
            iterator = tqdm(iterator)
        for i in iterator:
            mt   = self._mt_list[i]
            soff = self._stape_off[i]
            slen = self._stape_len[i]
            for j in range(slen):
                data_found = False
                mf         = self._stape[j + soff]
                current_mt = mt
                while current_mt is not None:
                    okey = 1000 * current_mt + mf
                    if okey in offset:
                        data_found = True
                        break
                    else:
                        current_mt = MTHierarchy[current_mt].parent()
                assert data_found, error('Transition data not found for MT={}, MF={}'.format(mt, mf))
                cc_pos = offset_keys.index(okey)
                stape  = mf + (cc_pos << 5)  # pack 5 bit MF + 27 bit position
                self._stape[j + soff] = stape

    def _prepareUnifiedMultiplicityAndDeposition(self, verbose: bool):
        if verbose:
            print(info('Generate global multiplicity & deposition table ...'))
        ngn = self._desc.ngn()
        egn = self._desc.egn()
        egg = self._desc.egg()
        self._multiplicity = np.empty(len(self._mt_list) * ngn, dtype=float)
        self._edepo        = np.empty(len(self._mt_list) * ngn, dtype=float)
        egn_mean = (egn[1:] + egn[:-1]) * 0.5
        egg_mean = (egg[1:] + egg[:-1]) * 0.5
        for i in range(len(self._mt_list)):
            mseg   = np.zeros(ngn, dtype=float)
            eseg   = np.zeros(ngn, dtype=float)
            mt     = self._mt_list[i]
            soff   = self._stape_off[i]
            slen   = self._stape_len[i]
            stape  = self._stape[soff:soff + slen] 
            mf     = stape & 31
            mf_unique, mf_counts = np.unique(mf, return_counts=True)
            # hadron multiplicity
            mseg += np.sum(mf_counts[mf_unique != 16])

            # gamma multiplicity
            if 16 in mf_unique:
                current_mt  = mt
                gamma_found = False
                while current_mt is not None:
                    if current_mt not in self._reaction.keys():
                        current_mt = MTHierarchy[current_mt].parent()
                        continue
                    if 16 in self._reaction[current_mt].keys():
                        gamma_found  = True
                        reaction     = self._reaction[current_mt]
                        mseg += reaction[16].multiplicity()
                        break
                    else:
                        current_mt = MTHierarchy[current_mt].parent()
                assert gamma_found, error('Gamma data not found for MT={}'.format(mt))

            # deposition
            reaction = self._reaction[mt]
            if 26 in reaction.keys():  # Tabulated depo data exist
                if verbose:
                    print(info('Energy deposition data found for MT={}, {} reaction'.format(mt, REACTION_TYPE[mt])))
                for group in range(ngn):
                    mask        = reaction[26].probabilityMask(group)
                    eseg[group] = np.sum(egn_mean * mask)
                pass
            else:  # depo data not exist
                za  = self._desc.za()
                za += 1  # neutron absorption
                for mf_this, count_this in zip(mf_unique, mf_counts):
                    if mf_this in GENDF_MF_TO_ZA.keys():
                        za -= GENDF_MF_TO_ZA[mf_this] * count_this
                if not za:
                    if verbose:
                        print(info('Residual target remnant not found for MT={}, {} reaction. Energy deposition is set to 0.'.format(mt, REACTION_TYPE[mt])))
                    eseg[:] = 0.0
                else:
                    if verbose:
                        print(info('Energy deposition data not found for MT={}, {} reaction. Using Q-value to calculate deposition.'.format(mt, REACTION_TYPE[mt])))
                    q_value = reaction.qvalue()
                    eseg[:] = q_value
                    for group in range(ngn):
                        for mf_this, count_this in zip(mf_unique, mf_counts):
                            reaction_parent = reaction
                            if mf_this not in reaction.keys():
                                current_mt = mt
                                data_found = False
                                while current_mt is not None:
                                    if current_mt not in self._reaction.keys():
                                        current_mt = MTHierarchy[current_mt].parent()
                                        continue
                                    if mf_this in self._reaction[current_mt].keys():
                                        data_found = True
                                        break
                                    else:
                                        current_mt = MTHierarchy[current_mt].parent()
                                assert data_found, 'MF={} data not found for MT={}'.format(mf, mt)
                                reaction_parent = self._reaction[current_mt]
                            mask = reaction_parent[mf_this].probabilityMask(group)
                            if mf_this == 16:
                                eseg[group] -= np.sum(egg_mean * mask) * reaction_parent[mf_this].multiplicity()[group]
                            else:
                                eseg[group] -= np.sum(egn_mean * mask) * count_this
                    eseg = egn_mean + eseg
            eseg[eseg < 0.0] = 0.0
            self._multiplicity[i * ngn:(i + 1) * ngn] = mseg
            self._edepo[i * ngn:(i + 1) * ngn]        = eseg * 1e-6  # eV to MeV
               
    def _generateGroupAlias(self, verbose: bool):
        if verbose:
            print(info('Generate group alias table ...'))
        iterator = self._gcontrol
        if verbose:
            iterator = tqdm(iterator)
        for epos, _, elen in iterator:
            if elen < 0:
                continue
            seg   = self._gprob[epos: epos + elen]
            table = AliasTable(np.arange(elen), seg)
            self._galias[epos: epos + elen] = table.alias()
            self._gprob[epos: epos + elen]  = table.prob()

    def _generateEquiprobAngles(self, endf: ENDF, nebins: int, verbose: bool):
        if verbose:
            print(info('Generate equiprobable angular distribution ...'))
        legpoly_list = self._eabin
        self._eabin  = np.empty((self._eabin.shape[0], nebins + 1), dtype=float)
        modifier     = (np.arange(0, legpoly_list.shape[1], 1) * 2 + 1) / 2   # ENDF legendre coeff

        # prepare MT=2 and high energy range (fast elastic)
        # first reaction is always MT=2 elastic
        cpos_min, _, len_min = self._gcontrol[self._len_thermal]

        cpos_max, _, len_max = self._gcontrol[self._desc.ngn() - 1]
        cpos_max = cpos_max + len_max + 1

        egn = GENDF.__egn

        # fast elastic scattering
        awr   = self._desc_origin.mass()
        awr   = max(1.0, awr)
        alpha = (awr - 1)**2 / (awr + 1)**2

        ad_mt2 = endf[2].ad
        is_cm  = ad_mt2.isOnCMSystem()
        # generate elastic equiprob
        if verbose:
            print(info('Angular distribution from the elastic scattering'))
        iterator = range(self._len_thermal, self._desc.ngn())
        if verbose:
            iterator = tqdm(iterator)
        for gin in iterator:
            cpos, gmin, len_c = self._gcontrol[gin]

            if len_c < 0:
                continue

            # get transition prob
            pseg   = self._gprob[cpos: cpos + len_c]

            ein = (egn[gin] + egn[gin + 1]) * 0.5

            dist = ad_mt2.getDistribution(ein)

            # get mu bins from transition probability
            mu_min = -1.0
            mu_max = +1.0
            if not is_cm:
                mu_min = muCMToLab(awr, mu_min)
                mu_max = muCMToLab(awr, mu_max)

            ebin_seg_list, _ = dist.getEquibin(nebins, mu_min, mu_max, pseg)

            if is_cm:
                ebin_seg_list = muCMToLab(awr, ebin_seg_list)

            self._eabin[cpos:cpos + len_c] = ebin_seg_list

        if verbose:
            print(info('Angular distribution from the GENDF'))
        iterator = range(len(legpoly_list))
        if verbose:
            iterator = tqdm(iterator)
        for i in iterator:
            if cpos_min <= i < cpos_max:  # skip MT=2 high energy scattering
                continue
            poly           = legpoly_list[i]
            ebin_seg, _    = legendreToEquibin(poly * modifier, nebins)
            self._eabin[i] = ebin_seg
        return
    
    def write(self, file_name: str):
        file = Fortran(file_name, mode='w')
        # header
        file.write(np.array((self._desc.za()), dtype=np.int32))                     # ZA
        file.write(np.array((self._desc_origin.isomericNumber()), dtype=np.int32))  # isomeric number
        file.write(np.array((self._temperature), dtype=np.float32))                 # temperature
        file.write(np.fromstring(self._sab, dtype=np.uint8))                        # SAB name
        file.write(np.array((self._eabin.shape[1]), dtype=np.int32))                # number of angle group
        # XS
        file.write(self._reaction[1].xs().astype(np.float32))
        # reaction
        file.write(self._mt_list.astype(np.int32))
        file.write(self._ralias.flatten().astype(np.int32))
        file.write(self._rprob.flatten().astype(np.float32))
        file.write(self._stape_len.astype(np.int32))
        file.write(self._stape_off.astype(np.int32))
        file.write(self._stape.astype(np.uint32))
        # group
        file.write(self._multiplicity.astype(np.float32))
        file.write(self._edepo.astype(np.float32))
        file.write(self._gcontrol.flatten().astype(np.int32))
        file.write(self._galias.astype(np.int32))
        file.write(self._gprob.astype(np.float32))
        file.write(self._eabin.flatten().astype(np.float32))
        file.close()

    def keys(self):
        return self._reaction.keys()
    
    def __getitem__(self, index):
        return self._reaction[index]
    
