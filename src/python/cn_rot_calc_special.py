import numpy as np


def read_inp(inp_file):
    inp = open(inpfile, 'r')
    params = {}
    lines = []
    for line in inp:
            lines.append(line)
            params["nfile"]  = str(lines[1]).strip()
            params["ntimes"] = int(lines[3].split()[0])
            params["dt"]     =float(lines[3].split()[1])
            params["volume"] = float(lines[5])
            params["mol_name"] = str(lines[7]).strip()
            inp.close()

            params["L"] = params["volume"]**(1/3.)
    print("Input parameters read in")
    for key in params:
        print("%8s = %8s" % (key, str(params[key])))
    return params

def read_traj(f,t)
