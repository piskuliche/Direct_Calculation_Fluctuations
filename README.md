# Direct\_Calculation\_Fluctuations

This code acts as an automation script (and analysis script) for MD simulations run using the Direct Calculation Method

The github repo should be cloned into the directory where you want to run the simulation in.

From this point forward the following nomenclature will be used:

sim/ - top level directory, where simulation runs
sim/dcm - Git Repo Level
sim/dcm/src/ - source directory in the dcm folder



#sim/dcm/backbone.sh

This is the main code you interact with-
run by: bash backbone.sh

Required Input: input\_file (in sim/ directory)
Output: Quite a bit, these are just the ones that it specifically creates when you run it.

sim/dcm/.flag\_instruct
sim/in.nve
sim/in.water
sim/dcm/.flag\_compile
sim/dcm/exec/msd\_rot\_calc
sim/dcm/mol\_names
sim/dcm/.flag\_traj
sim/dcm/.flag\_innve
sim/dcm/.flag\_filesystem
sim/file\_names
sim/file\_setup.py
sim/set\_msd\_calcs.py
sim/nve.sh
sim/job\_array.sh
sim/sub.sh
sim/msd\_rot\_calc
sim/dcm/.flag\_nve
sim/dcm/.flag\_checkup
sim/logs/array\*
sim/logs/direct\_calc\_nve\*
sim/dcm/.flag\_grabflucts
sim/grab\_flucts.log
sim/dcm/.flag\_valcalc
sim/flucts\_calc.py
sim/test.inp
sim/logs/msd\_calc.o\*

#sim/dcm/src/input/read\_input.sh

Called: Everytime the backbone script is called

This script reads input file parameters into the backbone script. An example input\_file is copied to the sim/ directory in the first step of the backbone.sh script

#sim/dcm/src/python/file\_setup.py

required input: file\_names (created automatically by backbone.sh)
        mol\_names (created automatically by backbone.sh)
produced output: setup\_files (which produces a host of other files)
         sub\_script

setup\_files output:
sim/FILES
for i subdirectories
sim/FILES/$i
sim/FILES/$i/nve.sh
sim/FILES/$i/set\_msd\_calcs.py
sim/FILES/$i/msd\_rot\_calc
sim/FILES/$i/restart.$i

sub\_script output:
sim/job\_array\_[0-$j]


#sim/dcm/src/python/set\_msd\_calcs.py

legacy - includes a set number of times and timestep value

This code automatically creates the msd input file.
Run automatically with syntax: python set\_msd\_calc.py -inp # -mol molname

required input: log.lammps (produced by nve simulations)
produced output: msd\_rot\_calc.in (in each FILES subfolder)

#sim/dcm/src/python/checkup.py

This program just checks the log file of every sub simulation a
nd makes sure that it exists - if it doesnt then it outputs the commands ne
cessary to rerun the file.

required input: none
produced output: commands to rerun missing dirs (to screen)


#sim/dcm/src/python/flucts\_calc.py

lecacy = includes a set number of times (in msd arrays)
future changes - add column number to calculation input file read

This program does the actual calculation of the weighted correlation functions and produces the output of the correlation functions.

Run as: python flucts\_calc.py -inp input\_file\_name -files number of files -blocks #num blocks -mol molecule name

required input: input file (with names of energies)
                val\_init.out (for each value in input file)
produced output: tons, listed now for each value
    ea\_msd\_val\_molname.dat 
    d\_val\_molname\_msd\_tot.dat
    d2\_val\_molname\_msd\_tot.dat
    dc2\_val\_molname.dat
    d2c2\_val\_molname.dat
    c2\_total\_molname\_result.dat
    and for each block:
        bl\_\#\_molname\_c2\_val.dat
        bl\_\#\_molname\_dc2\_val\_val.dat
        bl\_val\_molname\_msd.dat
        bl\_val\_molname\_c2.dat

#sim/dcm/src/python/grab\_flucts.py

legacy = Need to add an input file (the same as flucts calc) that just provides the item and the location in the log file.

This is the program that actually looks through each of the log files and grabs each of the energies/volumes/whatever thermodynamic data you want to weight the correlation function by.

required input: log.lammps (of each nve trajectory)
produced output: val\_init.out (needed by flucts\_calc.py, for each energy/vol/etc)

#sim/dcm/src/sub/job\_array.sh

This is the job array that actually runs the nve runs for all directories and then calculates the msds for each one.

#sim/dcm/src/sub/nve.sh

This is the submission script that is used by checkup.py to resubmit jobs that haven't run for one reason or another.

#sim/dcm/src/sub/sub.sh

This is the script that subits the setup\_files file which builds the filesystem. This was created by file\_setup.py

#sim/dcm/src/dependencies

This includes an example lammps input file for the npt and nve runs involved in this, an example data file, and a torque submission script

#sim/dcm/src/fortran/msd\_rot\_calc.f90

legacy: this has a set number of dimensions of arrays at compilation, these will need to be resized eventually to be effective.
        There is also currently no setting that lets you move beyond 3 atom molecules - this should be fixed soon as well to be as open as possible. i.e. there needs to be an option to set the atoms\_per\_mol, and then choose how to calculate things.


This is thie msd calculation code that makes the magic happen. It calculates the msds and the reorientation correlation functions (c1 and c2).

required input: msd\_rot\_calc.in (created by set\_msd\_calc.py)
                traj\_#\_molname.xyz (created by lammps)
produced output: c1\_\#\_molname.dat
                 c2\_\#\_molname.dat
                 msd\_\#\_molname.dat
