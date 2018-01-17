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
    if [ -f ../file_names ]; then
        rm ../file_names
    fi
    for ((i=$startconfig; i<=$endconfig; i+=$sepconfig )); do 
        echo $i >> ../file_names
    done
    mkdir ../FILES
    cp src/python/file_setup.py ../
    cp src/python/set_msd_calcs.py ../
    cp src/sub/nve.sh ../
    cp src/sub/job_array.sh ../
    cp src/sub/sub.sh ../
    cp src/exec/msd_rot_calc ../
    cd ../
    python file_setup.py
    msub sub.sh
    cd -
    echo "  Filesystem is being built"
    echo "  Please Wait until The Job Completes Before Moving Forward"
    touch .flag_filesystem
    exit 1
fi

# STEP: RUN NVE TRAJECTORIES

FILE=.flag_nve
if [ -f $FILE ]; then
        echo "-NVE Flag Exists"
else
    echo "-NVE Flag Missing"
    echo "  Running NVE Trajectories"
    cd ../
    bash sub_script
    cd -
    echo "  NVE Trajectories Submitted"
    touch .flag_nve
    exit 1
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
    mkdir logs
    mv array* logs
    mv direct_calc_nve* logs
    cd -
    echo "  Wait until remaining NVE trajectories run"
    echo "  If this step doesn't work the first time,"
    echo "  check your input files"
    touch .flag_checkup
    exit 1
fi

#STEP: Grab Flucts
FILE=.flag_grabflucts
if [ -f $FILE ]; then
    echo "-Grab Flucts Flag Exists"
else
    echo "-Grab Flucts Flag Missing"
    echo "  Running Flucts Grab Code"
    cp src/python/grab_flucts.py ../
    cd ../
    python grab_flucts.py > grab_flucts.log
    cd -
    echo "  Grab Flucts Has Finished!"
    touch .flag_grabflucts
    exit 1
fi


# STEP: CALCULATE WEIGHTED CORRELATION FUNCTIONS/Act. Eners./Act. Vols
FILE=.flag_valcalc
if [ -f $FILE ]; then
    echo "-Fluctuation Calculation Flag Exist"
else
    echo "-Fluctuation Calculation Flag Missing"
    echo "  Running Fluctuation Calc"
    cp src/python/flucts_calc.py ../
    cp src/input/test.inp ../
    cd ../
    mv msd_calc.o* logs
    python flucts_calc.py -inp test.inp -files 5000 -blocks 10
    cd -
    echo "  Fluctuation Calculation Completed"
    touch .flag_valcalc
    exit 1
fi
