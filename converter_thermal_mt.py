
import os
import subprocess
import time

import numpy as np

from src.endf.endf import ENDF
from src.prompt import Prompt

def printHelp():
    print("Parameters: --input        | -i  <path>      Path of ENDF database                 ")
    print("            --scattering   | -s  <path>      Path of ENDF thermal scattering data  ")
    print("            --output       | -o  <path>      Path for RT2NDL results               ")
    print("            --threads      | -th <int>       Number of threads                     ")
    print("            --target       | -t  <path>      Thermal target file                   ")
    print("            --equiprobable | -e  <int>       Number of equiprobable angle bins     ")
    print("            --help         | -h              Print this message                    ")


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
    print("ENDF thermal scattering path must be specified")
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

argv = Prompt()["--threads", "-th"]
nthreads = 1
if argv:
    nthreads = int(argv[0])

argv = Prompt()["--target", "-t"]
target_path = ""
if not argv:
    print("Library out path must be specified")
    printHelp()
    exit(1)
else:
    target_path = argv[0]

argv = Prompt()["--equiprobable",  "-e"]
if not argv:
    print("Number of equiprobable angle bins must be specified")
    printHelp()
    exit(1)
nebins = int(argv[0])

# get ENDF header lists
input_list  = {}
scatt_list  = {}

# build input/output lists
for file in os.listdir(input_path):
    endf_data = ENDF(os.path.join(input_path, file), verbose=False, read_only_header=True)
    endf_desc = endf_data.desc()
    if endf_desc.isomericNumber():  # isomeric nucleus -> pass
        continue
    za = endf_desc.za()
    input_list[za] = file

for file in os.listdir(scatt_path):
    endf_data = ENDF(os.path.join(scatt_path, file), verbose=False, read_only_header=True)
    endf_desc = endf_data.desc()
    mat       = endf_desc.mat()
    sab       = endf_desc.zsymam().__repr__().strip().replace(' ', '')
    scatt_list[mat] = (file, sab)


# initialize subprocess
procs = [
    subprocess.Popen(['python', 'subprocess_empty.py'], stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
] * nthreads
current_target = [None] * nthreads

os.makedirs(output_path, exist_ok=True)

# prepare folders
for i in range(nthreads):
    os.makedirs('thread{}'.format(i), exist_ok=True)

# prepare log
logs = [None] * nthreads
for i in range(nthreads):
    logs[i] = open('thread{}.log'.format(i), mode='w')

# prepare target
target_list = []
with open(target_path) as file:
    for line in file:
        items = line.split()
        if len(items) == 0:
            continue
        if items[0][0] == '#':
            continue
        mat     = int(items[0])
        za      = int(items[1])
        temp    = float(items[2])
        mfactor = float(items[3])
        target_list += [[mat, za, temp, mfactor]]

while True:
    stalled = False
    for i, proc in enumerate(procs):
        if proc.poll() is None:
            stalled = True
            continue
        else:
            if proc.poll() != 0:
                print('Fail to convert ENDF MAT {} for isotope {} at {} K'.format(current_target[i][0], current_target[i][1], current_target[i][2]))
            if len(target_list) == 0:  #finished
                continue
            target = target_list.pop()
            current_target[i] = target
            input_file = os.path.join(input_path,  input_list[target[1]])
            scatt_file = os.path.join(scatt_path,  scatt_list[target[0]][0])
            # output syntax
            out_name   = "{}_{}K_{}.bin".format(target[1], int(round(target[2])), scatt_list[target[0]][1])
            out_file   = os.path.join(output_path, out_name)
            command    = [
                'python', 'converter_thermal.py', 
                '-i', input_file, 
                '-s', scatt_file,
                '-o', out_file, 
                '-w', 'thread{}'.format(i), 
                '-e', str(nebins), 
                '-t', str(target[2]),
                '-m', str(target[3])
            ]
            stalled = True
            print('Convert ENDF MAT {} for isotope {} at {} K'.format(target[0], target[1], target[2]))
            procs[i] = subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=logs[i], stderr=logs[i])
    if not stalled:
        break
    time.sleep(0.5)

for i in range(nthreads):
    logs[i].close()

exit(0)