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
    echo "-Please create an input_file, one is included in src/input"
    echo "-Please create an flucts.inp, one is included in src/input"
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
    pgf90 -O3 flux_side.f90 -o flux_side
    mv msd_rot_calc ../exec/
    mv visc_calc ../exec/
    mv flux_side ../exec/
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

# STEP: RUN NVE TRAJECTORIES

FILE=.flag_nve
if [ -f $FILE ]; then
        echo "-NVE Flag Exists"
else
    echo "-NVE Flag Missing"
    read -r -p "Would you like to setup and run the NVE trajectories? [y/N] " response
    if [[ "$response" =~ ^([yY][eE][sS]|[yY])+$ ]]
    then
        echo "  Setting up and running NVE Trajectories"
        if [ -f ../file_names ]; then
            rm ../file_names
        fi
        for ((i=$start_config; i<=$end_config; i+=$sep_config )); do
            echo $i >> ../file_names
        done
        mkdir ../FILES
        # Copies various codes up to main directory
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
        cp src/exec/flux_side ../
        cd ../
        # Generate Submission Script
        python gen_sub_scripts.py
        if [ $timestep = 'FALSE' ]; then
            python non-unif-sample.py
        else
            python unif-sample.py
        fi
        screen -d -m bash run_array 
        cd -  
        echo "  NVE Trajectories Submitted"
        touch .flag_nve
    else
        echo "  NVE Trajectories were NOT submitted."
    fi
    exit 0
fi

# Check for nve completion flag:
FILE=../.flag_nvecomplete
if [ -f $FILE ]; then
    echo "-NVE Completion Flag Exists"
else
    echo "-NVE Completion Flag Missing"
    echo " Please be patient and wait for NVE trajectories to finish."
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
            if [ $program = 'LAMMPS' ]; then
                tail -n1 $i/log.lammps | grep -q Total || echo "cd $i; msub nve.sh; cd ../" >> rerun
            elif [ $program = 'CP2K' ]; then
                echo "  Since you are using CP2K, you need to check individually whether runs have finished."
            fi
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
    msub grabfluctsub.sh
    cd -
    echo "  Grab Flucts is running as a job"
    echo "  Please be patient - and check that grab_flucts.py is completed before continuing"
    echo "  Grab Flucts Has Finished!"
    touch .flag_grabflucts
    exit 0
fi


# STEP: INITIALIZE SEGMENTS CORRELATION FUNCTIONS/Act. Eners./Act. Vols
FILE=.flag_initarray
if [ -f $FILE ]; then
    echo "-InitArray Calculation Flag Exist"
else
    echo "-InitArray Calculation Flag Missing"
    echo "  Running InitArray Calc"
    cp src/python/init_flucts.py ../
    cp src/python/do_flucts.py ../
    cp src/shell/init_array.sh ../
    cd ../
    mkdir SEG
    msub init_array.sh
    cd -
    echo "  Fluctuation Array Submitted"
    echo "  Please wait until all jobs are completed."
    touch .flag_initarray
    exit 0
fi


#STEP: Fluctuations Calculation
FILE=.flag_valcalc
if [ -f $FILE ]; then
    echo "-Fluctuations Calculation Flag Exists"
else
    echo "-Fluctuations Flag Missing"
    echo "  Running Fluctuations Calculation"
    cp src/shell/fluctssub.sh ../
    cd ../
        for ((x=1; x<=$num_molecs; x++))
        {
            if [ $x -eq 1 ]
            then
                echo "python do_flucts.py flucts.inp msd $molec1 $blocks" > fluctssub.sh
                echo "python do_flucts.py flucts.inp c2 $molec1 $blocks" >> fluctssub.sh
                if [ $program = 'LAMMPS' ]; then
                    echo "python do_flucts.py flucts.inp shear $molec1 $blocks" >> fluctssub.sh
                fi
            fi
            if [ $x -eq 2 ]
            then
                echo "python do_flucts.py flucts.inp msd $molec2 $blocks" >> fluctssub.sh
                echo "python do_flucts.py flucts.inp c2 $molec2 $blocks" >> fluctssub.sh
                if [ $program = 'LAMMPS' ]; then
                    echo "python do_flucts.py flucts.inp shear $molec2 $blocks" >> fluctssub.sh
                fi 
            fi
            if [ $x -eq 3 ]
            then
                echo "python do_flucts.py flucts.inp msd $molec3 $blocks" >> fluctssub.sh
                echo "python do_flucts.py flucts.inp c2 $molec3 $blocks" >> fluctssub.sh
                if [ $program = 'LAMMPS' ]; then
                    echo "python do_flucts.py flucts.inp shear $molec3 $blocks" >> fluctssub.sh
                fi 
            fi
            if [ $x -eq 4 ]
            then
                echo "python do_flucts.py flucts.inp msd $molec4 $blocks" >> fluctssub.sh
                echo "python do_flucts.py flucts.inp c2 $molec4 $blocks" >> fluctssub.sh
                if [ $program = 'LAMMPS' ]; then
                    echo "python do_flucts.py flucts.inp shear $molec4 $blocks" >> fluctssub.sh
                fi 
            fi
        }
    msub fluctssub.sh 
    cd -
    echo "  Fluctuations has finished!"
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
