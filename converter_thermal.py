
import os
import pathlib
import shutil
import subprocess
import time

import numpy as np

from src.prompt import Prompt
from src.njoy import NJOYInput
from src.settings import getSetting, ENV
from src.fortran import Fortran
from src.gendf.gendf import GENDF
from src.endf.endf import ENDF


def printHelp():
    print("Parameters: --input        | -i  <filename>    ENDF input                         ")
    print("            --scattering   | -s  <filename>    ENDF S(a,b) table                  ")
    print("            --output       | -o  <filename>    Output of Group-wised library      ")
    print("            --workspace    | -w  <path>        NJOY working directory             ")
    print("            --temperature  | -t  <float>       Temperature for Doppler broadening ")
    print("            --molecule     | -m  <float>       Molecular number factor            ")
    print("            --equiprobable | -e  <int>         Number of equiprobable angle bins  ")
    print("            --verbose      | -v                Activate verbose log               ")
    print("            --help         | -h                Print this message                 ")


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

argv = Prompt()["--scattering", "-s"]
scatt_path = ""
if not argv:
    print("S(a,b) input path must be specified")
    printHelp()
    exit(1)
else:
    scatt_path = argv[0]

argv = Prompt()["--output", "-o"]
output_path = ""
if not argv:
    print("Library out path must be specified")
    printHelp()
    exit(1)
else:
    output_path = argv[0]

argv = Prompt()["--workspace", "-w"]
working_directory = None
if argv:
    working_directory = argv[0]

argv = Prompt()["--temperature", "-t"]
temperature = 293.6
if argv:
    temperature = float(argv[0])

argv = Prompt()["--molecule", "-m"]
mfactor = 1.0
if argv:
    mfactor = float(argv[0])

argv = Prompt()["--equiprobable",  "-e"]
if not argv:
    print("Number of equiprobable angle bins must be specified")
    printHelp()
    exit(1)
nebins = int(argv[0])

argv    = Prompt()["--verbose", "-v"]
verbose = False
if argv is not None:
    verbose = True

getSetting("settings/setting.txt", ENV)

egn = Fortran(ENV["njoy_ngroup"]).read(np.float32)
egg = Fortran(ENV["njoy_ggroup"]).read(np.float32)

# MeV to eV
egn = egn * 1e6
egg = egg * 1e6

endf_data = ENDF(input_path, verbose=False, read_only_header=False)
endf_desc = endf_data.desc()

ad = endf_data[2].ad.getDistribution(1.0)

mat = endf_desc.mat()
za  = endf_desc.za()

scatt_data = ENDF(scatt_path, verbose=False, read_only_header=True)
scatt_desc = scatt_data.desc()

mat_thermal = scatt_desc.mat()
za_thermal  = scatt_desc.za()
sab_str     = scatt_desc.zsymam().__repr__()
sab_str     = sab_str.strip().replace(' ', '')

# detect
print("*** ENDF MAT={}, ISOTOPE ZA={} IS DETECTED ***".format(mat, za))

# njoy data processing
print("*** NJOY processing ***")
if working_directory:
    os.makedirs(working_directory, exist_ok=True)
else:
    working_directory = '.'

njoy_input_file  = ENV["njoy_input"]
njoy_output_file = ENV["njoy_output"]

njoy_target      = int(ENV["njoy_target"])
njoy_kernel      = int(ENV["njoy_kernel"])
njoy_result      = int(ENV["njoy_GENDF"])
njoy_target_file = os.path.join('tape{}'.format(njoy_target))
njoy_kernel_file = os.path.join('tape{}'.format(njoy_kernel))
njoy_result_file = os.path.join('tape{}'.format(njoy_result))
shutil.copy(input_path, os.path.join(working_directory, njoy_target_file))
shutil.copy(scatt_path, os.path.join(working_directory, njoy_kernel_file))

os.chdir(working_directory)
pathlib.Path(njoy_output_file).unlink(missing_ok=True)

ninput = NJOYInput(njoy_input_file)
ninput.setEnv(mat, temperature)
ninput.setNeutronGroup(egn)
ninput.setPhotonGroup(egg)
ninput.moder(njoy_target, -21)
ninput.reconr(-21, -22, 0.0005)
ninput.broadr(-21, -22, -23, 0.0005)
ninput.thermr(njoy_kernel, -23, -24, mat_thermal, 2, 1, 0.005, 4)
ninput.groupr(-21, -24, -30, 1, 1, 7, 8, 1e7)
ninput.moder(-30, njoy_result)
ninput.stop()
ninput.write()

with subprocess.Popen([ENV["njoy_executable"],
                       "-i", njoy_input_file,
                       "-o", njoy_output_file]) as process:
    while True:
        if process.poll() == 0:
            break
        elif process.poll() is None:
            time.sleep(0.2)
        else:
            raise OSError

os.chdir('..')

print("*** GENDF data processing ***")
gendf_data = GENDF(os.path.join(working_directory, njoy_result_file), nebins, temperature, endf_data, sab=sab_str, mfactor=mfactor, verbose=verbose)
print('*** Write data "{}" ***'.format(output_path))
gendf_data.write(output_path)

exit(0)