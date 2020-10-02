#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
from subprocess import PIPE, Popen

def solve_it(data):
    process = Popen(['./TSPSolver.exe'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    (stdout, stderr) = process.communicate(bytes(data, encoding="UTF-8"))
    if stderr is not None and len(stderr) != 0:
        print(stderr)
    return str(stdout, encoding="UTF-8")

import sys

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        file_location = sys.argv[1].strip()
        with open(file_location, 'r') as input_data_file:
            input_data = input_data_file.read()
        print(solve_it(input_data))
    else:
        print('This test requires an input file.  Please select one from the data directory. (i.e. python solver.py ./data/tsp_51_1)')