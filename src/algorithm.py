from __future__ import annotations

"""
Provides interpolation algorithm that follows ENDF interpolation law,
and equi-probable angular bin conversion algorithm from legendre polynomial

this code is part of the RT2 project
"""

import numpy as np
from scipy.special import legendre
from scipy.interpolate import interp1d, interp2d
from scipy.integrate import quad, trapezoid
from scipy.optimize import root_scalar

from src import constants as const


def legendreToPolynomial(coeff: np.ndarray) -> np.ndarray:
    """
    Convert legendre coefficient array to polynomial

    :param coeff: Legendre coefficients 
    :return: Polynomial coefficients
    """
    polynomial = np.zeros(len(coeff), dtype=np.float64)
    for order, val in enumerate(coeff):
        polynomial[-1-order:] += val * legendre(order)
    return polynomial


def legendreToEquibin(coeff: np.ndarray, nbin: int, mu_min: float=-1.0, mu_max: float=1.0) -> tuple:
    """
    Convert legendre angle distribution to
    equiprobable angle bin distribution
    from mu_min to mu_max

    :param coeff: legendre coefficient series
    :param nbin: the number of equiprobable angle bin
    :param mu_min: lower bound of directional cosine
    :param mu_max: upper bound of directional cosine
    :return: (equiprobable bin, total area of legendre polynomial)
    """
    if mu_min > mu_max:
        raise ValueError("mu_min must be smaller than mu_max")
        
    poly = legendreToPolynomial(coeff)

    # find roots, only real number
    roots = np.roots(poly)
    roots = np.real(roots[np.isreal(roots)])
    roots = roots[(mu_min < roots) * (roots < mu_max)]
    roots = np.sort(roots)
    roots = np.unique(roots)
    if mu_min not in roots:
        roots = np.append(mu_min, roots)
    if mu_max not in roots:
        roots = np.append(roots, mu_max)

    # get integral
    polyint = np.poly1d(np.polyint(poly))
    # get area between each neighboring roots
    area_cumul = polyint(roots)
    area       = area_cumul[1:] - area_cumul[:-1]
    area_total = np.sum(area[area > 0])
    area_seg   = area_total / nbin
    last_area  = 0

    # find equiprob angle bin
    angle_bin     = np.empty(nbin+1, dtype=np.float64)
    angle_bin[0]  = mu_min
    angle_bin[-1] = mu_max
    n = 1
    for i in range(len(area)):
        if area[i] <= 0:
            continue
        root_lower = roots[i]
        root_upper = roots[i+1]
        int_lower  = polyint(root_lower)
        while True: # get answer
            polyint_t = np.copy(polyint)
            polyint_t[-1] += -int_lower + last_area - n * area_seg
            roots_int = np.roots(polyint_t)
            roots_int = np.real(roots_int[np.isreal(roots_int)])
            roots_int = roots_int[(root_lower <= roots_int) * (roots_int <= root_upper)]
            if len(roots_int) > 0:
                angle_bin[n] = np.min(roots_int)
                n += 1
            else:
                break
        last_area += area[i]

    return angle_bin, area_total


def logMean(a: float | np.ndarray, b: float | np.ndarray) -> float | np.ndarray:
    """
    Get logarithm mean

    :param a: min value
    :param b: max value
    :return: logarithm mean of a and b
    """
    return (b - a) / (np.log(b) - np.log(a))



class interp1d:
    def __init__(self, x: np.ndarray, y: np.ndarray, inte: int):
        """
        one-dimensional interpolation, follows ENDF interpolation law

        x: domain of the function.
        y: abscissa of the function.
        inte: ENDF interpolation law.
        """
        self._int = inte
        if self._int == 2: # linear-linear
            self._f = interp1d(x, y)
        elif self._int == 3: # linear-log
            self._f = interp1d(x, np.log(y))
        elif self._int == 4: # log-linear
            self._f = interp1d(np.log(x), y)
        elif self._int == 5: # log-log
            self._f = interp1d(np.log(x), np.log_y)
        else:
            raise ValueError('illegal interpolation law')
    
    def get(self, x: float) -> float:
        """
        Get abscissa value at position x
        """
        if self._int == 2:
            y = self._f(x)
        elif self._int == 3:
            log_y = self._f(x)
            y = np.exp(log_y)
        elif self._int == 4:
            y = self._f(np.log(x))
        elif self._int == 5:
            log_y = self._f(np.log(x))
            y = np.exp(log_y)
        return y


class interp2d:
    def __init__(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, inte: int):
        """
        two-dimensional interpolation, follows ENDF interpolation law

        :param x: domain of the function.
        :param y: domain of the function.
        :param z: value of the function. Shape must be (nx, ny)
        :param inte: ENDF interpolation law.
        """
        self._int = inte
        if self._int == 2: # linear-linear
            self._f = interp2d(x, y, z)
        elif self._int == 3: # linear-log
            # 0 point handling
            z_copy = np.copy(z)
            z_copy[z_copy < const.LOG_MIN] = const.LOG_MIN
            self._f = interp2d(x, y, np.log(z_copy))
        elif self._int == 4: # log-linear
            x_copy = np.copy(x)
            y_copy = np.copy(y)
            x_copy[x_copy < const.LOG_MIN] = const.LOG_MIN
            y_copy[y_copy < const.LOG_MIN] = const.LOG_MIN
            self._f = interp2d(np.log(x_copy), np.log(y_copy), z)
        elif self._int == 5: # log-log
            z_copy = np.copy(z)
            z_copy[z_copy < const.LOG_MIN] = const.LOG_MIN
            x_copy = np.copy(x)
            y_copy = np.copy(y)
            x_copy[x_copy < const.LOG_MIN] = const.LOG_MIN
            y_copy[y_copy < const.LOG_MIN] = const.LOG_MIN
            self._f = interp2d(np.log(x_copy), np.log(y_copy), np.log(z_copy))
        else:
            raise ValueError('illegal interpolation law')

    def get(self, x: float, y: float) -> float:
        """
        Get value at (x, y)
        """
        if self._int == 2:
            z = self._f(x, y)
        elif self._int == 3:
            log_z = self._f(x, y)
            z = np.exp(log_z)
        elif self._int == 4:
            z = self._f(np.log(x), np.log(y))
        elif self._int == 5:
            log_z = self._f(np.log(x), np.log(y))
            z = np.exp(log_z)
        return z


def getInterpFtnCumulArea(xx: np.ndarray, yy: np.ndarray, x: float) -> float:
    """
    get the area of numerical function (xx, yy) from xx[0] to x

    :param xx: domain of the function.
    :param yy: abscissa of the function.
    :param x: upper bound of domain.
    :return: total area from xx[0] to x
    """
    if x < xx[0] or x > xx[-1]:
        raise ValueError("x must be in xx range")
    target_index = np.argmax(x <= xx)
    y = interp1d(xx, yy, 2).get(x)
    xx_new = np.append(xx[:target_index], x)
    yy_new = np.append(yy[:target_index], y)
    return trapezoid(yy_new, xx_new)


def getInterpFtnCumulValue(xx: np.ndarray, yy: np.ndarray, area: float) -> float:
    """
    get the upper bound x when the integrated area of
    numerical function (xx, yy) from xx[0] to x is the 'area'

    :param xx: domain of the function
    :param yy: value of the function
    :param area: area: area of integral
    :return: upper bound x
    """    
    interpolation = interp1d(xx, yy, kind='linear', fill_value="extrapolate")
    
    def integral_to_x(upper_x):
        result, _ = quad(interpolation, xx[0], upper_x)
        return result

    def objective(upper_x):
        return integral_to_x(upper_x) - area

    solution = root_scalar(objective, bracket=[xx[0], xx[-1]], method='bisect')
    
    if solution.converged:
        return solution.root
    else:
        raise ValueError('Cannot find solution x')



class AliasTable:
    def __init__(self, domain, prob):
        """
        convert discrete probability density function to
        alias table

        :param domain: domain of the probability density function
        :param prob: value of the probability density function
        """
        self._domain = domain
        self._alias_table = -np.ones(domain.shape, dtype=int)
        self._prob_table  = np.copy(prob)
        prob_tag = np.ones(prob.shape, dtype=bool)

        assert len(prob), 'Alias table got empty probablity table'

        mean = np.sum(prob) / len(prob)
        # set alias table
        for i in range(len(prob) - 1):
            lower = np.where((self._prob_table < mean) * prob_tag)[0]
            upper = np.where((self._prob_table > mean) * prob_tag)[0]

            if len(lower) == 0 or len(upper) == 0:
                continue

            target_low = lower[0]
            target_up = upper[0]

            aux = mean - self._prob_table[target_low]
            self._prob_table[target_up]  -= aux
            self._prob_table[target_low] /= mean
            self._alias_table[target_low] = target_up

            prob_tag[target_low] = False

        self._prob_table[prob_tag] = 10

    def sampling(self):
        rand = np.random.random()
        aj = rand * len(self._domain)
        j = int(aj)
        aj -= j
        if aj > self._prob_table[j]:
            j = self._alias_table[j]
        return j

    def alias(self):
        return self._alias_table

    def prob(self):
        return self._prob_table

def probFromAlias(alias_table: np.ndarray, prob_table: np.ndarray) -> np.ndarray:
    """
    convert alias table to the probability density function
    
    :param alias_table: Alias table
    :param prob_table: probablility table

    :return: probability density function
    """
    prob = np.zeros(prob_table.shape)
    mean = 1 / len(prob_table)
    for i in range(len(prob_table)):
        if prob_table[i] >= 1:
            prob[i] += mean
        else:
            target        = alias_table[i]
            prob[target] += mean * (1.e0 - prob_table[i])
            prob[i]      += mean * prob_table[i]
    
    return prob
