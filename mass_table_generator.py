
import os
import pathlib
import shutil
import subprocess
import time

from tqdm import tqdm
import numpy as np
from pyne.endf import Evaluation

from src.prompt import Prompt
from src.fortran import Fortran
from src.constants import MASS_NEUTRON_DA

def printHelp():
    print("Parameters: --input        | -i  <path>      Path of ENDF database               ")
    print("            --output       | -o  <filename>  File name of mass table             ")
    print("            --help         | -h              Print this message                  ")


# Prompt
argv = Prompt()["--help", "-h"]
if argv is not None:
    printHelp()
    exit(1)

argv = Prompt()["--input", "-i"]
input_path = ""
if not argv:
    print("ENDF input path must be specified")
    printHelp()
    exit(1)
else:
    input_path = argv[0]

argv = Prompt()["--output", "-o"]
output_path = ""
if not argv:
    print("Library out path must be specified")
    printHelp()
    exit(1)
else:
    output_path = argv[0]


DTYPE_ATOM_LIST = [
    ('za',   np.int32  ), 
    ('mass', np.float32)
]


target_list = os.listdir(input_path)
atom_list   = {}
for i, file in enumerate(target_list):
    file      = os.path.join(input_path, file)
    endf_data = Evaluation(file, verbose=False)
    za        = endf_data.target['ZA']
    mass      = endf_data.target['mass'] * MASS_NEUTRON_DA
    atom_list[za] = mass

mass_array = np.empty(len(atom_list.keys()), DTYPE_ATOM_LIST)
for i, key in enumerate(sorted(atom_list)):
    mass_array[i]['za']   = key
    mass_array[i]['mass'] = atom_list[key]

mass_file = Fortran(output_path, mode='w')
mass_file.write(mass_array['za'])
mass_file.write(mass_array['mass'])
mass_file.close()

exit(0)