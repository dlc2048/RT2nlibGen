from __future__ import annotations

"""
Write NJOY input file

this code is part of the RT2 project
"""

import numpy as np


def _joinCard(*args):
    return " " + " ".join(list(map(str, args))) + " /\n"


def _joinGroup(arr):
    out = ""
    i = 0
    for data in arr:
        out += " "
        out += "{:6e}".format(data)
        if i == 3:
            out += "\n"
            i = 0
            continue
        i += 1
    out += " /\n"
    return out


class NJOYInput:
    def __init__(self, file_name):
        """
        NJOY input file
        visit NJOY21 homepage to get more information
        [https://www.njoy21.io/]
        """
        self._file = open(file_name, mode="w")
        self.mat = None
        self.temperature     = None
        self._egn = None
        self._egg = None

    def setEnv(self, mat, temperature):
        self.mat = mat
        self.temperature = temperature

    def setNeutronGroup(self, egn: np.ndarray):
        self._egn = egn
    
    def setPhotonGroup(self, egg: np.ndarray):
        self._egg = egg

    def moder(self, nin, nout):
        """writle MODER module input"""
        # header
        self._file.write("moder\n")
        # cards
        self._file.write(_joinCard(nin, nout))

    def reconr(self, nin, nout, err):
        """writle RECONR module input"""
        # header
        self._file.write("reconr\n")
        # cards
        self._file.write(_joinCard(nin, nout))
        self._file.write(_joinCard("'pendf tape for" + str(self.mat) + "'"))
        self._file.write(_joinCard(self.mat))
        self._file.write(_joinCard(err))
        self._file.write(_joinCard(0))

    def broadr(self, nendf, nin, nout, err):
        """writle BROADR module input"""
        # header
        self._file.write("broadr\n")
        # cards
        self._file.write(_joinCard(nendf, nin, nout))
        self._file.write(_joinCard(self.mat))
        self._file.write(_joinCard(err))
        self._file.write(_joinCard(self.temperature))
        self._file.write(_joinCard(0))

    def thermr(self, nendf, nin, nout, kernel_mat, iin, icoh, tol, emax):
        """writle THERMR module input"""
        # header
        self._file.write("thermr\n")
        # cards
        self._file.write(_joinCard(nendf, nin, nout))
        self._file.write(_joinCard(kernel_mat, self.mat, 10, 1,
                                   iin, icoh, 0, 1, 221, 1))
        self._file.write(_joinCard(self.temperature))
        self._file.write(_joinCard(tol, emax))

    def groupr(self, nendf, npend, ngout2, ign, igg, iwt, lord, sigz, mt_lists: None | list = None):
        """writle GROUPR module input"""
        # header
        self._file.write("groupr\n")
        # cards
        self._file.write(_joinCard(nendf, npend, 0, ngout2))
        self._file.write(_joinCard(self.mat, ign, igg, iwt, lord))
        self._file.write(_joinCard("'neutron group structure of " + str(self.mat) + "'"))
        self._file.write(_joinCard(self.temperature))
        self._file.write(_joinCard(sigz))
        if self._egn is not None:
            self._file.write(_joinCard(len(self._egn) - 1))
            self._file.write(_joinGroup(self._egn))
        if self._egg is not None:
            self._file.write(_joinCard(len(self._egg) - 1))
            self._file.write(_joinGroup(self._egg))
        # target reactions
        self._file.write(_joinCard(3))
        self._file.write(_joinCard(3, 221))
        self._file.write(_joinCard(6))
        self._file.write(_joinCard(6, 221))

        for i in range(21, 26):
            self._file.write(_joinCard(i))
        if mt_lists is None:
            self._file.write(_joinCard(26))  # all
        elif 19 not in mt_lists:
            self._file.write(_joinCard(26))
        self._file.write(_joinCard(12))
        self._file.write(_joinCard(16))
        for i in range(2):
            self._file.write(_joinCard(0))

    def moder(self, nin, nout):
        """writle MODER module input"""
        # header
        self._file.write("moder\n")
        # cards
        self._file.write(_joinCard(nin, nout))

    def stop(self):
        """writle STOP module input"""
        # header
        self._file.write("stop")
        
    def write(self):
        self._file.close()


