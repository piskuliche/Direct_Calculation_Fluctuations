"""
This is the program to create the:
    "input_file"
    "flucts.inp"
"""

def int_input(paramname,desc):
    try:
        print("%s" % desc)
        param = int(raw_input('Please enter an integer for %s\n' % paramname))
    except ValueError:
        print("Entry is not an integer")
        sys.exit(1)
    return param

def float_input(paramname,desc):
    try: 
        print("%s" % desc)
        param = float(raw_input('Please enter a float for %s\n' % paramname))
    except ValueError:
        print("Entry is not a float")
        sys.exit(1)
    return param

def str_input(paramname,desc):
    try:
        print("%s" % desc)
        param = str(raw_input('Please enter a string for %s\n' % paramname))
    except ValueError:
        print("Entry is not a string")
        sys.exit(1)
    return param

def par_input(type, paramname,desc):
    if type == "int":
        param = int_input(paramname,desc)
        return param
    elif type == "float":
        param = float_input(paramname,desc)
        return param
    elif type == "str":
        param = str_input(paramname,desc)
        return param
    else:
        print("Incorrect type error - Cannot generate input")
        sys.exit(1)


input = 'input_file'
inp = open(input, 'w')

#define input types
typei = "int"
typef = "float"
types = "str"

inp.write("#------------------------------------\n")
inp.write("#- This is the input file for:      -\n")
inp.write("#-   Direct Fluctions Calculations  -\n")
inp.write("#- Only change uncommented lines    -\n")
inp.write("#------------------------------------\n")

desc = "Program (LAMMPS or CP2K):"
inp.write("# %s \n" % desc)
paramname = "progname"
progname = par_input(types,paramname,desc)
inp.write("%s\n" % progname)

desc = "Run Style (NVT or NPT):"
inp.write("# %s\n" % desc)
paramname = "runstyle"
runstyle=par_input(types, paramname, desc)
inp.write("%s\n" % runstyle)

desc = "Number of Molecule Types:"
inp.write("# %s\n" % desc)
paramname = "nummoltype"
nummoltype = par_input(typei, paramname, desc)
inp.write("%s\n" % nummoltype)

desc = "Identifier of Molecules(LEAVE 4 LINES SPACE):"
inp.write("# %s\n" % desc)
for i in range(4):
    if i <= nummoltype-1:
        desc = "Molecule name"
        molname = par_input(types, paramname, desc)
        inp.write("%s\n" % molname)
    else:
        inp.write("\n")

desc = "Starting Configuration (in STEPS):"
inp.write("# %s\n" % desc)
paramname = "startconfig"
startconfig = par_input(typei,paramname,desc)
inp.write("%s\n" % startconfig)

desc = "Ending Configuration (in STEPS):"
inp.write("# %s\n" % desc)
paramname = "endconfig"
endconfig = par_input(typei,paramname,desc)
inp.write("%s \n" % endconfig)

desc = "Separation of Configurations (in STEPS):"
inp.write("# %s\n" % desc)
paramname = "sepconfig"
sepconfig = par_input(typei, paramname, desc)
inp.write("%s \n" % sepconfig)

desc = "# Dump Frequency (in picoseconds):"
inp.write("# %s \n" % desc)
paramname = "dumpfreq"
dumpfreq = par_input(typef, paramname, desc)
inp.write("%s \n" % dumpfreq)

desc = "# Num Times (How many dumps there are per nve trajectory - (NVE Length/dump freq):"
inp.write("# %s\n" % desc)
paramname = "dumpnum"
dumpnum = par_input(typei, paramname, desc)
inp.write("%s \n" % dumpnum)

desc = "Blocks (How many blocks there should be):"
inp.write("# %s\n" % desc)
paramname = "numblocks"
numblocks = par_input(typei, paramname, desc)
inp.write("%s \n" % numblocks)

desc = "Files (How many NVE trajectories there are):"
inp.write("# %s\n" % desc)
paramname = "numfiles"
numfiles = par_input(typei, paramname, desc)
inp.write("%s\n" % numfiles)

desc = "NVE Length (Number of STEPS):"
inp.write("# %s\n" % desc)
paramname = "nvesteps"
nvesteps = par_input(typei, paramname, desc)
inp.write("%s\n" % nvesteps)

desc = "Correlation Functions (TRANSPORT or IONPAIRING):"
inp.write("# %s\n" % desc)
paramname = "corrfunc"
corrfunc = par_input(types, paramname, desc)
inp.write("%s\n" % corrfunc)
