from __future__ import annotations

import re

"""
Interpret and parse ENDF document
also provides particle id and reaction type comments

this code is part of the RT2 project
"""

MF_TYPE = {
    1 : "General Information"             ,
    3 : "Reaction Cross Sections"         ,
    4 : "Angular Distributions"           ,
    5 : "Energy Distributions"            ,
    6 : "Energy-Angle Distributions"      ,
    12: "Photon Production Yield Data"    , 
    13: "Photon Production Cross Sections"
}


#### reaction type ####
REACTION_TYPE = {
    1  : "(total)"     , 2  : "(elastic)"   , 3  : "(nonelastic)",
    4  : "(z,n')"      , 5  : "(any)"       , 10 : "(continuum)" ,
    11 : "(z,2nd)"     , 16 : "(z,2n)"      , 17 : "(z,3n)"      ,
    18 : "(fission)"   , 19 : "(z,f)"       , 20 : "(z,nf)"      ,
    21 : "(z,2nf)"     , 22 : "(z,na)"      , 23 : "(z,n3a)"     ,
    24 : "(z,2na)"     , 25 : "(z,3na)"     , 27 : "(absorp)"    ,
    28 : "(z,np)"      , 29 : "(z,n2a)"     , 30 : "(z,2n2a)"    ,
    32 : "(z,nd)"      , 33 : "(z,nt)"      , 34 : "(z,nHe-3)"   ,
    35 : "(z,nd2a)"    , 36 : "(z,nt2a)"    , 37 : "(z,4n)"      ,
    38 : "(z,3nf)"     , 41 : "(z,2np)"     , 42 : "(z,3np)"     ,
    44 : "(z,n2p)"     , 45 : "(z,npa)"     , 91 : "(z,nc)"      ,
    101: "(disapp)"    , 102: "(z,r)"       , 103: "(z,p)"       ,
    104: "(z,d)"       , 105: "(z,t)"       , 106: "(z,He-3)"    ,
    107: "(z,a)"       , 108: "(z,2a)"      , 109: "(z,3a)"      ,
    111: "(z,2p)"      , 112: "(z,pa)"      , 113: "(z,t2a)"     ,
    114: "(z,d2a)"     , 115: "(z,pd)"      , 116: "(z,pt)"      ,
    117: "(z,da)"      , 201: "(z,xn)"      , 202: "(z,xr)"      ,
    203: "(z,xp)"      , 204: "(z,xd)"      , 205: "(z,xt)"      ,
    206: "(z,xHe-3)"   , 207: "(z,xa)"      , 221: "(thermal)"   , 
    649: "(z,pc)"      , 699: "(z,dc)"      , 749: "(z,tc)"      , 
    799: "(z,3-Hec)"   , 849: "(z,ac)"      , 891: "(z,2nc)"     
}

REACTION_TYPE[451] = "(z,...)"

for i in range(41): # (z,n') reactions
    REACTION_TYPE[i+50] = "(z,n" + str(i) + ")"

for i in range(49): # (z,p') reactions
    REACTION_TYPE[i+600] = "(z,p" + str(i) + ")"

for i in range(49): # (z,d') reactions
    REACTION_TYPE[i+650] = "(z,d" + str(i) + ")"

for i in range(49): # (z,t') reactions
    REACTION_TYPE[i+700] = "(z,t" + str(i) + ")"

for i in range(49): # (z,He-3') reactions
    REACTION_TYPE[i+750] = "(z,3-He" + str(i) + ")"

for i in range(49): # (z,a') reactions
    REACTION_TYPE[i+800] = "(z,a" + str(i) + ")"

for i in range(26): # (z,2n') reactions
    REACTION_TYPE[i+875] = "(z,2n" + str(i) + ")"


MT_REPR_PATTERN = r"\(z,(?:(\d*)n[^,]*)?(?:(\d*)p[^,]*)?(?:(\d*)d[^,]*)?(?:(\d*)t[^,]*)?(?:(\d*)He-3[^,]*)?(?:(\d*)a[^,]*)?\)"

# Multiplicity
NEUTRON_MULT  = {}
PROTON_MULT   = {}
DEUTERON_MULT = {}
TRITIUM_MULT  = {}
HE3_MULT      = {}
ALPHA_MULT    = {}

for key in REACTION_TYPE:
    value = REACTION_TYPE[key]
    match = re.fullmatch(MT_REPR_PATTERN, value)
    if match:
        n   = int(match.group(1)) if match.group(1) else (1 if 'n'    in value else 0)
        p   = int(match.group(2)) if match.group(2) else (1 if 'p'    in value else 0)
        d   = int(match.group(3)) if match.group(3) else (1 if 'd'    in value else 0)
        t   = int(match.group(4)) if match.group(4) else (1 if 't'    in value else 0)
        he3 = int(match.group(5)) if match.group(5) else (1 if '3-He' in value else 0)
        a   = int(match.group(6)) if match.group(6) else (1 if 'a'    in value else 0)
    else:
        n   = 0
        p   = 0
        d   = 0
        t   = 0
        he3 = 0
        a   = 0
    NEUTRON_MULT[key]  = n
    PROTON_MULT[key]   = p
    DEUTERON_MULT[key] = d
    TRITIUM_MULT[key]  = t
    HE3_MULT[key]      = he3
    ALPHA_MULT[key]    = a

NEUTRON_MULT[2]  = 1  # elastic
NEUTRON_MULT[18] = 1  # fission

# MT hierarchy (See sum rules for cross sections in Section 0.5.1 Table 14).
class MTHierarchyNode:
    def __init__(self, mt: int):
        self._mt     = mt
        self._parent = None
        self._child  = []

    def mt(self) -> int:
        return self._mt

    def parent(self) -> int | None:
        return self._parent
    
    def child(self) -> list:
        return self._child
    
    def setParent(self, parent: int):
        self._parent = parent

    def addChild(self, child: int):
        if child not in self._child:
            self._child += [child]


class MTHierarchy:
    __node = {}

    @staticmethod
    def init():
        MTHierarchy.__node = {}
        for mt in REACTION_TYPE.keys():
            MTHierarchy.__node[mt] = MTHierarchyNode(mt)

    @staticmethod
    def addNodeLink(parent: int, child: int):
        assert(parent != child)
        parent_node = MTHierarchy.__node[parent]
        child_node  = MTHierarchy.__node[child]

        if (child_node.parent() is not None) and (child_node.parent() != parent):
            assert(False)
        
        child_node.setParent(parent)
        parent_node.addChild(child)

    @staticmethod
    def nodes():
        return MTHierarchy.__node
    
    @staticmethod
    def __class_getitem__(index: int):
        return MTHierarchy.__node[index]
    

MTHierarchy.init()

# Total
MTHierarchy.addNodeLink(1, 2)
MTHierarchy.addNodeLink(1, 3)

# Non-elastic
__NON_ELASTIC = [4, 5, 11, 16, 17, 18]
for mt in range(22, 26):
    __NON_ELASTIC += [mt]
for mt in range(28, 31):
    __NON_ELASTIC += [mt]
for mt in range(33, 38):
    __NON_ELASTIC += [mt]
__NON_ELASTIC += [41, 42, 44, 45]
for mt in __NON_ELASTIC:
    MTHierarchy.addNodeLink(3, mt)

# Neutron inelastic
for mt in range(50, 92):
    MTHierarchy.addNodeLink(4, mt) 

# (n,2n)
for mt in range(875, 892):
    MTHierarchy.addNodeLink(16, mt)

# fission
for mt in (19, 20, 21, 38):
    MTHierarchy.addNodeLink(18, mt)  

# absorption
MTHierarchy.addNodeLink(27, 101)

# disappearance
for mt in range(102, 110):
    MTHierarchy.addNodeLink(101, mt)
for mt in range(111, 118):
    MTHierarchy.addNodeLink(101, mt)


# (n,p)
for mt in range(600, 650):
    MTHierarchy.addNodeLink(103, mt)

# (n,d)
for mt in range(650, 700):
    MTHierarchy.addNodeLink(104, mt)

# (n,t)
for mt in range(700, 750):
    MTHierarchy.addNodeLink(105, mt)

# (n,He3)
for mt in range(750, 800):
    MTHierarchy.addNodeLink(106, mt)

# (n, alpha)
for mt in range(800, 850):
    MTHierarchy.addNodeLink(107, mt)
