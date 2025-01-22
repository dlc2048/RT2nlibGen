from __future__ import annotations

"""
GPUMC-CNDL global settings

this code is part of the RT2 project
"""


ENV = {}


def getSetting(file_name, setting_dict):
    with open(file_name) as file:
        lines = file.readlines()
        file.close()
    for line in lines:
        if len(line) == 0:
            continue
        if line[0] == "*": # comment line
            continue
        item = line.split()
        if len(item) == 0:
            continue
        setting_dict[item[0]] = item[1]
