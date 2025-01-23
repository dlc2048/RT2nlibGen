
import os
import pathlib
import shutil
import subprocess
import time

import numpy as np
from tqdm import tqdm
from pyne.endf import Evaluation

from src.prompt import Prompt

def printHelp():
    print("Parameters: --input        | -i  <path>      Path of ENDF database                 ")
    print("            --output       | -o  <path>      Path for RT2NDL results               ")
    print("            --threads      | -th <int>       Number of threads                     ")
    print("            --temperature  | -t  <float>     Temperature for Doppler broadening [K]")
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

argv = Prompt()["--temperature", "-t"]
temperature = 293.6
custom_temp = False
if argv:
    temperature = float(argv[0])
    custom_temp = True

argv = Prompt()["--equiprobable",  "-e"]
if not argv:
    print("Number of equiprobable angle bins must be specified")
    printHelp()
    exit(1)
nebins = int(argv[0])


# build input/output lists
target_list = os.listdir(input_path)
output_list = []
for file in target_list:
    file      = os.path.join(input_path, file)
    endf_data = Evaluation(file, verbose=False)
    za        = endf_data.target['ZA']
    meta      = endf_data.target['isomeric_state']
    za_str    = str(za)
    if meta:
        za_str += 'm{}'.format(meta)
    oformat = '{}_{}.bin'.format(za_str, int(np.round(temperature)))
    output_list += [oformat]

# initialize subprocess
procs = [
    subprocess.Popen(['python', 'subprocess_empty.py'], stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
] * nthreads
current_target = [''] * nthreads

# prepare folders
for i in range(nthreads):
    os.makedirs('thread{}'.format(i), exist_ok=True)

while True:
    stalled = False
    for i, proc in enumerate(procs):
        if proc.poll() is None:
            stalled = True
            continue
        else:
            if proc.poll() != 0:
                print('Fail to convert ENDF file {}'.format(current_target[i]))
            if len(target_list) == 0:  #finished
                continue
            target = target_list.pop()
            out    = output_list.pop()
            current_target[i] = target
            target_file = os.path.join(input_path, target)
            out_file    = os.path.join(output_path, out)
            command     = [
                'python', 'converter.py', 
                '-i', target_file, 
                '-o', out_file, 
                '-w', 'thread{}'.format(i), 
                '-e', str(nebins), 
                '-t', str(temperature)
            ]
            stalled = True
            print('Convert ENDF {}'.format(target))
            procs[i] = subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if not stalled:
        break
    time.sleep(0.5)

exit(0)