#!/bin/bash

#This is the backbone program for the direct calculation method.
#Currently it can calculate:
#   -Activation Energies
#   -Activation Volumes
#of Diffusion and Reorientation. 

module load compiler/pgi/15

echo "Running direct calculations."

# STEP: Compile MSD_ROT_CALC
FILE=.flag_compile    
if [ -f $FILE ]; then
    echo "-Compilation Flag Exists"
else
    echo "-Compilation Flag Missing"
    echo "  Compiling MSD_ROT_CALC"
    cd src/fortran
    pgf90 -O3 msd_rot_calc.f90 -o msd_rot_calc
    mv msd_rot_calc ../exec/
    cd ../../
    echo "  Compilation complete"
    touch .flag_compile
fi

# STEP: RUN TRAJECTORY
FILE=.flag_traj
if [ -f $FILE ]; then
    echo "-Trajectory Run Flag Exists"
else
    echo "  Trajectory Run Flag Missing"
    if [ -f "../log.lammps" ]; then
        echo "Trajectory Found"
        touch .flag_traj
    else
        echo "  At this point, you should run the base trajectory"
        echo "  Example NVT and NPT input files are located in src/dependencies"
        touch ../log.lammps # delete later - for testing only
        exit 1
    fi
fi

# STEP: CREATE NVE INPUT FILE
FILE=.flag_innve
if [ -f $FILE ]; then
    echo "-NVE Input Flag Exists"
else
    echo "-NVE Input Flag Missing"
    if [ -f "../in.nve" ]; then
        echo "  Input File Found"
        touch .flag_innve
    else
        echo "  At this point, you should provide an nve simulation input file"
        echo "  An example NVE input file is located in src/dependencies"
        echo "  Place this file in the main run directory to continue"
        touch ../in.nve # delete later - for testing only
        exit 1
    fi
fi

# STEP: FILE SYSTEM SETUP
FILE=.flag_filesystem
if [ -f $FILE ]; then
    echo "-Filesystem Flag Exists"
else
    echo "-Building filesystem"
    echo "  Please enter starting configuration, then hit [ENTER]"
    read startconfig
    echo "  Please enter ending configuration, then hit [ENTER]"
    read endconfig
    echo "  Please enter separation of configurations, then hit [ENTER]"
    read sepconfig
    
    if [[ -z "$startconfig" ]]; then
        printf "No input entered"
        exit 1
    else
        # If userInput is not empty show what the user typed in and run ls -l
        echo "You entered" "$startconfig"
    fi

    if [[ -z "$endconfig" ]]; then
        echo "No input entered"
        exit 1
    else
        # If userInput is not empty show what the user typed in and run ls -l
        echo "You entered" "$endconfig"
    fi

    if [[ -z "$sepconfig" ]]; then
        echo "No input entered"
        exit 1
    else
        # If userInput is not empty show what the user typed in and run ls -l
        echo "You entered " "$sepconfig"
    fi

    for ((i=$startconfig; i<$endconfig; i+=$sepconfig )); do 
        echo $i >> ../file_names
        touch ../RESTART/restart.$i
    done
    mkdir ../FILES
    cp src/python/file_setup.py ../
    cp src/sub/water_nve.sh ../
    cp src/sub/job_array.sh ../
    cd ../
    python file_setup.py
    bash setup_files
    cd -
    echo "-Filesystem is built"
    touch .flag_filesystem
fi

# STEP: RUN NVE TRAJECTORIES

FILE=.flag_nve
if [ -f $FILE ]; then
        echo "-NVE Flag Exists"
else
    echo "-NVE Flag Missing"
    echo "  Running NVE Trajectories"
    bash sub_script
    echo "  NVE Trajectories Submitted"
    touch .flag_nve
fi

# STEP: CHECK THAT NVE TRAJECTORIES RAN

FILE=.flag_checkup
if [ -f $FILE ]; then
    echo "-NVE Trajectory Checkup Flag Exists"
else
    echo "-NVE Trajectory Checkup Flag Missing"
    echo "  Running checkup.py"
    cp src/python/checkup.py ../
    cd ../
    python checkup.py > out
    bash out
    cd -
    echo "  Wait until remaining NVE trajectories run"
    echo "  If this step doesn't work the first time,"
    echo "  check your input files"
    touch .flag_checkup
fi

# STEP: CALCULATE MSD, C1, C2

FILE=.flag_ac
if [ -f $FILE ]; then
    echo "-Correlation Calculation Flag Exists"
else
    echo "-Correlation Calculation Flag Missing"
    echo "  Running MSD calculation"
     
fi
# STEP: CALCULATE WEIGHTED CORRELATION FUNCTIONS/Act. Eners./Act. Vols
# STEP: TABULATE RESULTS
