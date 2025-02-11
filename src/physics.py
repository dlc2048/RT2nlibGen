from __future__ import annotations

import numpy as np
from scipy import integrate, optimize


def energyToMuCM(alpha: float, inc_energy: float, out_energy: float):
    return (2*out_energy/inc_energy - (1+alpha))/(1-alpha)


def muCMToLab(awr: float, mu: float | np.ndarray) -> float | np.ndarray:
    x = 1 + awr * mu
    y = np.sqrt(awr**2 + 1 + 2 * awr * mu)
    return np.divide(x, y, out=np.zeros_like(x), where=y!=0)


def muLabToCM(awr: float, mu: float | np.ndarray) -> float | np.ndarray:
    return (mu**2 - 1 + mu * np.sqrt(awr**2 + mu**2 - 1)) / awr