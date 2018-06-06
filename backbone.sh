#!/bin/bash

#This is the backbone program for the direct calculation method.
#Currently it can calculate:
#   -Activation Energies
#   -Activation Volumes
#of Diffusion and Reorientation

module load compiler/pgi/15

echo "Running direct calculations."

# STEP: Initial Instructions
FILE=.flag_instruct
if [ -f $FILE ]; then
    echo "-Instructions Flag Exists"
else
    echo "-Instructions Flag Doesn't Exist"
    echo "-Please modify input_file"
    cp src/input/input_file ../
    cp src/input/flucts.inp ../
    touch .flag_instruct
    exit 0
fi

cp src/shell/read_input.sh ./
source read_input.sh
rm read_input.sh
# STEP: Compile MSD_ROT_CALC
FILE=.flag_compile    
if [ -f $FILE ]; then
    echo "-Compilation Flag Exists"
else
    echo "-Compilation Flag Missing"
    echo "  Compiling MSD_ROT_CALC"
    cd src/fortran
    pgf90 -O3 msd_rot_calc.f90 -o msd_rot_calc
    pgf90 -O3 visc_calc.f90 -o visc_calc
    mv msd_rot_calc ../exec/
    mv visc_calc ../exec/
    cd ../../
    echo "  Compilation complete"
    for ((x=1; x<=$num_molecs; x++))
        {
            if [ $x -eq 1 ]
            then
                echo $molec1 > mol_names
            fi
            if [ $x -eq 2 ]
            then
                echo $molec2 >> mol_names
            fi
            if [ $x -eq 3 ]
            then
                echo $molec3 >> mol_names
            fi
            if [ $x -eq 4 ]
            then
                echo $molec4 >> mol_names
            fi   
        }
    mv mol_names ../
    touch .flag_compile
fi

# STEP: RUN TRAJECTORY
FILE=.flag_traj
if [ -f $FILE ]; then
    echo "-Trajectory Run Flag Exists"
else
    echo "  Trajectory Run Flag Missing"
    if [ $program = 'LAMMPS' ]; then
        if [ -f "../log.lammps" ]; then
            echo "Trajectory Found"
            touch .flag_traj
        else
            echo "  At this point, you should run the base trajectory"
            echo "  Example NVT and NPT input files are located in src/dependencies"
            exit 0
        fi
    elif [ $program = 'CP2K' ]; then
        if [ -f "../cp2k.log" ]; then
            echo "Trajectory Found"
            mkdir ../RESTART
            mv ../*.restart ../RESTART
            cd ../RESTART
            echo "  What is the identifier for your restart files?"
            echo "  i.e. anything before the # in the .restart files"
            read identifier
            if [ -f "$identifier"$start_config".restart" ]; then
                echo "renaming restart files"
                for (( i=$start_config; i<=$end_config; i+=$sep_config ))
                    {
                        mv "$identifier"$i".restart" restart.$i
                    }
                touch .flag_traj
                cd -
            else
                echo "  Improper identifier, try again."
                exit 0
            fi
        else
            echo "  At this point, you should run the base trajectory"
            echo "  Example NVT and NPT input files are located in src/dependencies"
            exit 0
        fi
    else
        echo "  Incorrect program type in input_file"
        exit 0
    fi
fi

# STEP: CREATE NVE INPUT FILE
FILE=.flag_innve
if [ -f $FILE ]; then
    echo "-NVE Input Flag Exists"
else
    echo "-NVE Input Flag Missing"
    if [ $program = 'LAMMPS' ]; then
        if [ -f "../in.nve" ]; then
            echo "  Input File Found"
            touch .flag_innve
        else
            echo "  At this point, you should provide an nve (title in.nve) simulation input file"
            echo "  An example NVE input file is located in src/dependencies"
            echo "  Place this file in the main run directory to continue"
            exit 0
        fi
    elif [ $program = 'CP2K' ]; then
        if [ -f "../in.nve.cp2k" ]; then
            echo "  Input File Found"
            touch .flag_innve
        else
            echo "  At this point, you should provide an nve (title in.nve.cp2k) simulation input file"
            echo "  An example NVE input file is located in src/dependencies"
            echo "  Place this file in the main run directory to continue"
            exit 0
        fi 
    else
        echo "  Program information incorrect, check input_file"
        exit 0
    fi
fi

# STEP: FILE SYSTEM SETUP
FILE=.flag_filesystem
if [ -f $FILE ]; then
    echo "-Filesystem Flag Exists"
else
    echo "-Building filesystem"
    if [ -f ../file_names ]; then
        rm ../file_names
    fi
    for ((i=$start_config; i<=$end_config; i+=$sep_config )); do 
        echo $i >> ../file_names
    done
    mkdir ../FILES
    # Copies various codes up to simulation directory
    cp src/python/file_setup.py ../ 
    cp src/python/read_input.py ../
    cp src/python/set_msd_calcs.py ../
    cp src/python/gen_sub_scripts.py ../
    cp src/python/grab_press.py ../
    if [ $timestep = 'FALSE' ]; then
        cp src/python/non-unif-sample.py ../
    else
        cp src/python/unif-sample.py ../
    fi
    cp src/exec/msd_rot_calc ../
    cp src/exec/visc_calc ../
    cd ../
    # Generate Submission Scripts
    python gen_sub_scripts.py
    if [ $timestep = 'FALSE' ]; then
        python non-unif-sample.py 
    else
        python unif-sample.py 
    fi
    # Find and Replace in job_array.sh
    python file_setup.py
    msub sub.sh
    cd -
    echo "  Filesystem is being built"
    echo "  Please Wait until The Job Completes Before Moving Forward"
    echo "  Please also modify the script job_array.sh and nve.sh for your system"
    touch .flag_filesystem
    exit 0
fi


# STEP: RUN NVE TRAJECTORIES

FILE=.flag_nve
if [ -f $FILE ]; then
        echo "-NVE Flag Exists"
else
    echo "-NVE Flag Missing"
    read -r -p "Would you like to run the NVE trajectories? [y/N] " response
    if [[ "$response" =~ ^([yY][eE][sS]|[yY])+$ ]]
    then
        echo "  Running NVE Trajectories"
        cd ../
        screen -d -m bash sub_script
        cd -
        echo "  NVE Trajectories Submitted"
        touch .flag_nve

    else
        echo "  NVE Trajectories were NOT submitted."
    fi
    exit 0
fi

# STEP: CHECK THAT NVE TRAJECTORIES RAN

FILE=.flag_checkup
if [ -f $FILE ]; then
    echo "-NVE Trajectory Checkup Flag Exists"
else
    echo "-NVE Trajectory Checkup Flag Missing"
    echo "  Running checkup"
    cd ../FILES/
    # Checks and runs missing trajectories
    if [ -f rerun ]; then
        rm rerun
    fi 
    for (( i=$start_config; i<=$end_config; i+=$sep_config ))
        {
            tail -n1 $i/log.lammps | grep -q Total || echo "cd $i; msub nve.sh; cd ../" >> rerun
        } 
    
    if (( $(wc rerun | awk '{print $1}') > 200 )); then 
        echo "  Error: More than 100 folders to rerun - check ../FILES/rerun before continuing"
    else
        echo "  Submitting $(wc rerun | awk '{print $1}') missed jobs"
        bash rerun
    fi
    # Does some cleanup commands
    cd ../
    mkdir logs
    mv array* logs
    mv direct_calc_nve* logs
    cd Direct_Calculation_Fluctuations

    echo "  Wait until remaining NVE trajectories run"
    echo "  If this step doesn't work the first time,"
    echo "  check your input files"
    touch .flag_checkup
    exit 0
fi

#STEP: Grab Flucts
FILE=.flag_grabflucts
if [ -f $FILE ]; then
    echo "-Grab Flucts Flag Exists"
else
    echo "-Grab Flucts Flag Missing"
    echo "  Running Flucts Grab Code"
    cp src/python/grab_flucts.py ../
    cp src/shell/grabfluctsub.sh ../
    cd ../
    rm grab_flucts.log
    msub grabfluctsub.sh
    cd -
    echo "  Grab Flucts is running as a job"
    echo "  Please be patient - and check that grab_flucts.py is completed before continuing"
    echo "  Grab Flucts Has Finished!"
    touch .flag_grabflucts
    exit 0
fi


# STEP: CALCULATE WEIGHTED CORRELATION FUNCTIONS/Act. Eners./Act. Vols
FILE=.flag_valcalc
if [ -f $FILE ]; then
    echo "-Fluctuation Calculation Flag Exist"
else
    echo "-Fluctuation Calculation Flag Missing"
    echo "  Running Fluctuation Calc"
    cp src/python/flucts_calc.py ../
    cp src/shell/fluctssub.sh ../
    cd ../
    mv msd_calc.o* logs
    for ((x=1; x<=$num_molecs; x++))
        {
            if [ $x -eq 1 ]
            then
                echo "python flucts_calc.py -inp flucts.inp -files $num_files -blocks $blocks -mol $molec1 -ntimes $num_times -ind \$MOAB_JOBARRAYINDEX >> array_\$MOAB_JOBARRAYINDEX.o" >> fluctssub.sh
            fi
            if [ $x -eq 2 ]
            then
                echo "python flucts_calc.py -inp flucts.inp -files $num_files -blocks $blocks -mol $molec2 -ntimes $num_times -ind \$MOAB_JOBARRAYINDEX >> array_\$MOAB_JOBARRAYINDEX.o" >> fluctssub.sh
            fi
            if [ $x -eq 3 ]
            then
                echo "python flucts_calc.py -inp flucts.inp -files $num_files -blocks $blocks -mol $molec3 -ntimes $num_times -ind \$MOAB_JOBARRAYINDEX >> array_\$MOAB_JOBARRAYINDEX.o" >> fluctssub.sh
            fi
            if [ $x -eq 4 ]
            then
                echo "python flucts_calc.py -inp flucts.inp -files $num_files -blocks $blocks -mol $molec4 -ntimes $num_times -ind \$MOAB_JOBARRAYINDEX >> array_\$MOAB_JOBARRAYINDEX.o" >> fluctssub.sh
            fi
        }
    msub fluctssub.sh 
    cd -
    echo "  Fluctuation Calculation Completed"
    touch .flag_valcalc
    exit 0
fi

#STEP: Cleanup
FILE=.flag_cleanup
if [ -f $FILE ]; then
    echo "-Cleanup Flag Exists"
    echo "-There are no further options"
else
    echo "-Cleanup Flag Missing"
    echo "  Running cleanup"
    cd ../
    mkdir Analysis
    mv *msd.dat Analysis
    mv *c1.dat Analysis
    mv *c2.dat Analysis
    mv *shear.dat Analysis
    cd -
    echo "  Cleanup Has Finished!"
    touch .flag_cleanup
    exit 0
fi
