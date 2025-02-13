from __future__ import annotations

from src.endf.desc import NEUTRON_MULT, PROTON_MULT, DEUTERON_MULT, TRITIUM_MULT, HE3_MULT, ALPHA_MULT


GENDF_MF_TYPE = {1  : "info"         , 3  : "cross-section", 6  : "neutron"     , 16 : "gamma"       , 
                 21 : "proton"       , 22 : "deuteron"     , 23 : "tritium"     , 24 : "helium-3"    , 
                 25 : "alpha"        , 26 : "heavy-nuc"    }
GENDF_MF_SYMBOL = {6: "n", 16: "g", 21: "p", 22: "d", 23: "t", 24: "He-3", 25: "a"}

GENDF_MF_TO_MULT = {6 : NEUTRON_MULT, 21: PROTON_MULT, 22: DEUTERON_MULT, 
                    23: TRITIUM_MULT, 24: HE3_MULT,    25: ALPHA_MULT   }
GENDF_MF_TO_ZA   = {6: 1, 21: 1001, 22: 1002, 23: 1003, 24: 2003, 25: 2004}
GENDF_ZA_TO_MF   = {1: 6, 1001: 21, 1002: 22, 1003: 23, 2003: 24, 2004: 25}