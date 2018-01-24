#!/bin/bash

#------------------------------------
#- This is the input file for:      -
#-   Direct Fluctions Calculations  -
#- Only change uncommented lines    -
#------------------------------------
# Program (Lammps or CP2K):
program=`sed -n '7p' input_file`
# Run Style (NVT or NPT):
run_style=`sed -n '9p' input_file`
# Number of Molecules
num_molecs=`sed -n '11p' input_file`
# Identifier of Molecules
molec1=`sed -n '13p' input_file`
molec2=`sed -n '14p' input_file`
molec3=`sed -n '15p' input_file`
molec4=`sed -n '16p' input_file`
# Starting Configuration (in steps)
start_config=`sed -n '18p' input_file`
# Ending Configuration (in steps)
end_config=`sed -n '20p' input_file`
# Separation of Configurations (in steps)
sep_config=`sed -n '22p' input_file`
# Timestep (in ps)
timestep=`sed -n '24p' input_file`
# Blocks
blocks=`sed -n '26p' input_file`

echo "Input File Parameters"
echo "Program: "$program
echo "Run Style: "$run_style
echo "Num Molecs: "$num_molecs
echo "mol1: "$molec1
echo "mol2: "$molec2
echo "mol3: "$molec3
echo "mol4: "$molec4
echo "Start Config: "$start_config
echo "End Config: "$end_config
echo "Sep Config: "$sep_config
echo "Timestep: "$timestep
echo "Blocks: "$blocks
echo "End Input File Parameters"
