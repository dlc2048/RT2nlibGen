from __future__ import annotations

"""
ENDF FILE 4: ANGULAR DISTRIBUTIONS OF SECONDARY PARTICLES

this code is part of the RT2 project
"""

import numpy as np
import scipy.integrate as spi
import scipy.optimize as spo

from src.endf.stream import FileInterface, ENDFIfstream
from src.endf.record import RecCont
from src.endf.param import ANGULAR_DIST_TYPE, INTERPOLATION_SCHEME
from src.algorithm import Interp1d, legendreToEquibin


class DistributionFunction:
    def __init__(self, ftn: np.polynomial.Legendre | tuple):
        self._mode = isinstance(ftn, np.polynomial.Legendre)

        self._legendre = None
        self._interp   = None

        if self._mode:
            self._legendre = ftn
        else:
            self._interp = Interp1d(ftn[0], ftn[1], inte=2)

    def isLegendre(self):
        return self._mode
    
    def _getEquibinFromLegendre(self, nbin: int, mu_min: float, mu_max: float, pseg: np.ndarray) -> tuple:
        if mu_min > mu_max:
            raise ValueError("mu_min must be smaller than mu_max")
        
        coeff = self._legendre.coef

        # drop small coeffs
        if coeff[0] > 0.0:
            coeff /= coeff[0] * 2  # normalize
        coeff[np.abs(coeff) < 1e-10] = 0

        ftn = np.polynomial.Legendre(coeff)

        # find roots, only real number
        roots = ftn.roots()
        roots = np.real(roots[np.isreal(roots)])
        roots = roots[(mu_min < roots) * (roots < mu_max)]
        roots = np.sort(roots)
        roots = np.unique(roots)
        if mu_min not in roots:
            roots = np.append(mu_min, roots)
        if mu_max not in roots:
            roots = np.append(roots, mu_max)

        # get integral
        ftn_integ = ftn.integ(1)
        # get area between each neighboring roots
        area_cumul = ftn_integ(roots)
        area       = area_cumul[1:] - area_cumul[:-1]

        # equiprob area
        area_total = np.sum(area[area > 0])
        area_seg   = pseg * area_total / nbin

        angle_bin = np.zeros((len(pseg), nbin + 1))
        angle_bin[0,0]   = roots[0]
        angle_bin[-1,-1] = roots[-1]

        pcumul   = np.cumsum(pseg)
        pcum_seg = pcumul * area_total
        pcum_seg = np.append(0.0, pcumul)
        last_area = 0
        t = 0
        n = 1
        for i in range(len(area)):
            if area[i] <= 0:
                continue
            root_lower = roots[i]
            root_upper = roots[i+1]
            int_lower  = ftn_integ(root_lower)
            while True:  # get answer
                ftn_integ_t = ftn_integ - int_lower + last_area - n * area_seg[t] - pcum_seg[t]

                roots_int = ftn_integ_t.roots()
                roots_int = np.real(roots_int[np.isreal(roots_int)])
                roots_int = roots_int[(root_lower <= roots_int) * (roots_int <= root_upper)]
                if len(roots_int) > 0:
                    angle_bin[t,n] = np.min(roots_int)
                    n += 1
                else:
                    break
                
                if n == nbin + 1:
                    # inherit last bin
                    n  = 1
                    t += 1
                    if t >= len(pseg):
                        break
                    angle_bin[t, 0] = angle_bin[t - 1, -1]

            last_area += area[i]
            
        return angle_bin, area_total
    
    def _getEquibinFromTabula(self, nbin: int, mu_min: float, mu_max: float, pseg: np.ndarray) -> tuple:
        if mu_min > mu_max:
            raise ValueError("mu_min must be smaller than mu_max")
        
        x = self._interp.x()
        y = self._interp.y()
        
        mu_min = max(mu_min, x[0])
        mu_max = min(mu_max, x[-1])

        n0 = np.argmax(mu_min <= x)
        n1 = np.argmax(mu_max <= x)

        xn = x[n0:n1]
        if xn[0] > mu_min:
            xn = np.append(mu_min, xn)
        if xn[-1] < mu_max:
            xn = np.append(xn, mu_max)
        
        yn = np.empty_like(xn)
        for i in range(len(xn)):
            yn[i] = self._interp.get(xn[i])
        
        area       = 0.5 * (xn[1:] - xn[:-1]) * (yn[1:] + yn[:-1])
        area_total = np.sum(area)
        area_cum   = np.cumsum(area)

        angle_bin = np.zeros((len(pseg), nbin + 1))
        angle_bin[0,0]   = mu_min
        angle_bin[-1,-1] = mu_max

        pcumul    = np.cumsum(pseg)
        pcum_seg  = pcumul * area_total
        pcum_seg  = np.append(0.0, pcumul)

        t = 0
        n = 1
        for ni, acum in enumerate(area_cum):
            target_area = (1 - n / nbin) * pcum_seg[t] + n / nbin * pcum_seg[t + 1]
            if target_area > acum:
                continue
            print(1)

        return

    def getEquibin(self, nbin: int, mu_min: float, mu_max: float, pseg: np.ndarray) -> tuple:
        if self._mode:
            return self._getEquibinFromLegendre(nbin, mu_min, mu_max, pseg)
        else:
            return self._getEquibinFromTabula(nbin, mu_min, mu_max, pseg)


class AngularDist(FileInterface):
    def __init__(self, head: RecCont, stream: ENDFIfstream, mt: int):
        super().__init__(4, mt)
        self._ltt = ANGULAR_DIST_TYPE(head.l2())
        cont = stream.cont()
        self._li  = bool(cont.l1())
        self._lct = cont.l2() == 2  # data are given in the CM system if true
        self._legendre  = None
        self._tabulated = None
        if self._ltt == ANGULAR_DIST_TYPE.ALL_ISOTROPIC:
            pass
        elif self._ltt == ANGULAR_DIST_TYPE.LEGENDRE:
            self._legendre   = stream.tab2(True)
        elif self._ltt == ANGULAR_DIST_TYPE.TABULATED:
            self._tabulated = stream.tab2(False)
        elif self._ltt == ANGULAR_DIST_TYPE.MIXED:
            self._legendre  = stream.tab2(True)
            self._tabulated = stream.tab2(False)
        self._checkSectionEnd(stream)

    def isOnCMSystem(self) -> bool:
        return self._lct

    def _getDistributionLegendre(self, inc_e: float) -> np.polynomial.Legendre:
        tab = self._legendre.tab()
        for i, eq in enumerate(tab):
            e = eq.c2()
            if inc_e < e:
                break
        i0  = i - 1
        i1  = i
        e0  = tab[i0].c2()
        e1  = tab[i1].c2()
        eq0_base = tab[i0].list()
        eq1_base = tab[i1].list()
        dim_max = max(len(eq0_base), len(eq1_base))
        eq0 = np.zeros(dim_max)
        eq1 = np.zeros(dim_max)
        eq0[:len(eq0_base)] = eq0_base
        eq1[:len(eq1_base)] = eq1_base
        # do interpolation
        for iint, int in self._legendre.interp():
            if i0 < iint:
                break
        int = INTERPOLATION_SCHEME(int)
        if int == INTERPOLATION_SCHEME.LINEAR_LINEAR or int == INTERPOLATION_SCHEME.LINEAR_LOG:  # linear
            pass
        else:  # log
            inc_e = np.log(inc_e)
            e0    = np.log(e0)
            e1    = np.log(e1)
        coeff = (inc_e - e0) / (e1 - e0) * eq0 + (e1 - inc_e) / (e1 - e0) * eq1
        # append a0
        coeff    = np.append(1.0, coeff)
        modifier = (np.arange(0, len(coeff), 1) * 2 + 1) / 2   # ENDF legendre coeff
        return DistributionFunction(np.polynomial.Legendre(coeff * modifier))
    
    def _getDistributionTabulated(self, inc_e: float) -> tuple:
        tab = self._tabulated.tab()
        for i, eq in enumerate(tab):
            e = eq.c2()
            if inc_e < e:
                break
        i0  = i - 1
        i1  = i
        e0  = tab[i0].c2()
        e1  = tab[i1].c2()
        eq0_base   = tab[i0].tab()
        eq1_base   = tab[i1].tab()
        # do interpolation
        for iint, inte in self._tabulated.interp():
            if i0 < iint:
                break
        inte = INTERPOLATION_SCHEME(inte)
        interp_eq0 = Interp1d(eq0_base[:,0], eq0_base[:,1], inte=2)
        interp_eq1 = Interp1d(eq1_base[:,0], eq1_base[:,1], inte=2)

        domain = np.unique(np.concatenate((eq0_base[:,0], eq1_base[:,0])))
        value  = np.zeros_like(domain)

        for i, x in enumerate(domain):
            xarr = np.array((e0, e1))
            yarr = np.array((interp_eq0.get(x), interp_eq1.get(x)))
            interp = Interp1d(xarr, yarr, inte.value)
            value[i] = interp.get(inc_e)

        return DistributionFunction((domain, value))
    
    def _getDistributionMixed(self, inc_e: float) -> np.polynomial.Legendre | tuple:
        tab = self._legendre.tab()
        for i, eq in enumerate(tab):
            e = eq.c2()
            if inc_e < e:
                return self._getDistributionLegendre(inc_e)
        return self._getDistributionTabulated(inc_e)

    def getDistribution(self, inc_e: float) -> np.polynomial.Legendre | tuple:
        if self._ltt == ANGULAR_DIST_TYPE.ALL_ISOTROPIC:  # isotropic
            return np.polynomial.Legendre([0.5])
        elif self._ltt == ANGULAR_DIST_TYPE.LEGENDRE:
            return self._getDistributionLegendre(inc_e)
        elif self._ltt == ANGULAR_DIST_TYPE.TABULATED:
            return self._getDistributionTabulated(inc_e)
        elif self._ltt == ANGULAR_DIST_TYPE.MIXED:
            return self._getDistributionMixed(inc_e)
        assert False
    