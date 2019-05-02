#!/usr/bin/env python

import numpy as np
import os,sys,subprocess
import src.python.path

def step_machinename():
    exists = os.path.isfile('machine.name')
    if exists:
        f = open('machine.name', 'r')
        mach_name = f.readline().strip()
        homepath = f.readline().strip()
        f.close()
    else:
        print("->machine.name is not set!")
        sys.exit()
    return mach_name, homepath

def step_instruct():
    flag='.flag_instruct'
    flag_exists = os.path.isfile(flag)
    if flag_exists:
        print("Instructions Flag exists.")
    else:
        print("Instructions Flag doesn't exist")
        print("->Please create input_file")
        print("->Please create flucts.inp")
        print("->Alternatively, you can run the input_generator to create input_file")
        inpresponse = str(input("Would you like to do so? y/n\n"))
        print(inpresponse)
        if 'y' in inpresponse:
            import src.python.gen_input
        elif 'Y' in inpresponse:
            import src.python.gen_input
        open(flag, 'a').close()
    return


def step_traj():
    flag='.flag_traj'
    flag_exists = os.path.isfile(flag)
    if flag_exists:
        print("Trajectory Run Flag Exists.")
    else:
        print("Trajectory Run Flag Missing.")
        f = open('mol_names','w')
        for item in inputparam.molec:
            f.write('%s\n' % item)
        f.close()
        if inputparam.prog == "LAMMPS":
            log_exists = os.path.isfile('log.lammps')
            if log_exists:
                print("->Trajectory found")
                open(flag, 'a').close()
            else: 
                print("-->At this point, you should run the base trajectory")
                print("-->Example NVT and NPT input files are located in src/dependencies")
                sys.exit()
        elif inputparam.prog == "CP2K":
            log_exists = os.path.isfile('cp2k.log')
            if log_exists:
                print("->Trajectory found")
                if not os.path.exists("RESTART"):
                    os.makedirs("RESTART")
                subprocess.call("mv *.restart RESTART", shell=True)
                print("->What is the identifier for your restart files?")
                identifier=str(input("-->i.e. anything before the # in the .restart files"))
                if os.path.isfile(identifier+str(inputparam.start_config)+".restart"):
                    print("->Renaming restart files")
                    for i in range(inputparam.start_config,inputparam.end_config,inputparam.sep_config):
                        os.rename('RESTART/'+identifier+str(i)+'.restart','RESTART/restart.'+str(i))
                    open(flag, 'a').close()
                else:
                    print("-->Improper identifier, try again.")
                    print("-->If you would like to skip, please type touch .flag_traj")
            else:
                print("->At this point, you should run the base trajectory")
                print("->Example NVT and NPT input files are located in src/dependencies")
                sys.exit()
        else:
            print("Incorrect program type")
            sys.exit()
    return 

def step_nve():
    """
    This step makes the file_names file, the FILES directory, and then calls a couple of python
    codes to create submission scripts etc
    """
    if not os.path.isfile('in.nve') and not os.path.isfile('in.nve.cp2k'):
        print("Need to put an in.nve or in.nve.cp2k file into directory")
        sys.exit()
    for i in range(inputparam.num_molecs):
        if inputparam.molec[i] != "water" and not os.path.isfile(str(inputparam.molec[i])+'.txt'):
            print("Need to create %s.txt" % inputparam.molec[i])
    flag='.flag_nve'
    flag_exists = os.path.isfile(flag)
    if flag_exists:
        print("NVE Flag Exists.")
    else:
        print("NVE Flag Missing.")
        print("->Setting up and running NVE trajectories")
        f = open("file_names", 'w')
        # Writes the filenames
        for i in range(inputparam.start_config, inputparam.end_config+inputparam.sep_config, inputparam.sep_config):
            f.write("%d\n" % i)
        f.close()
        if not os.path.exists("FILES"):
            os.makedirs("FILES")
        import src.python.gen_sub_scripts
        if inputparam.timestep == 'FALSE':
            import src.python.non_unif_sample
        else:
            import src.python.unif_sample
        print("->NVE Trajectories are ready to be submitted as sbatch run_array.sh")
        open(flag, 'a').close()
        sys.exit()
    return

def step_nvecomplete():
    """
    This step exists to make sure that you don't go to the next real step 
    if simulations are still running. It basically just makes you double check
    that all jobs have finished.
    """
    flag='.flag_nvecomplete'
    flag_exists = os.path.isfile(flag)
    if flag_exists:
        print("NVE Completion Flag Exists.")
    else:
        print("NVE Completion Flag Missing.")
        print("->Note: you have to create this flag yourself")
        print("->Check the queue - make sure ALL jobs have finished")
        print("->Then, once you are satisfied type 'touch .flag_nvecomplete'")
        sys.exit()
    return

def step_checkup():
    flag='.flag_checkup'
    flag_exists = os.path.isfile(flag)
    if flag_exists:
        print("NVE Trajectory Checkup Flag Exists.")
    else:
        print("NVE Trajectory Checkup Flag Missing.")
        print("-> Running Checkup")
        f = open('rerun', 'w')
        for i in range(inputparam.start_config, inputparam.end_config+inputparam.sep_config, inputparam.sep_config):
            if (i%1000000) == 0:
                print("-->%d has been checked" % i)
            if os.path.isfile("msd_%d_%s.dat" % (i,inputparam.molec[0])):
                print("%d does not exist" % i)
                if inputparam.prog == "LAMMPS":
                    f.write("mkdir FILES/%d; cp in.nve FILES/%d; cd FILES/%d; sed -i -e 's@/panfs/pfs.local/scratch/thompson/e924p726/edit_direct_calc/Direct_Calculation_Fluctuations@%d@g' nve.sh; sbatch nve.sh; cd ../../\n" % (i,i,i,i))
                elif inputparam.prog == "CP2K":
                    f.write("mkdir FILES/%d; cp in.nve.cp2k FILES/%d; cd FILES/%d; sed -i -e 's@/panfs/pfs.local/scratch/thompson/e924p726/edit_direct_calc/Direct_Calculation_Fluctuations@%d@g' nve.sh; sbatch nve.sh; cd ../../\n" % (i,i,i,i))
        f.close()
        if os.path.isfile('rerun') and os.path.getsize('rerun'):
            f.open('rerun','r')
            count = 0
            for line in f:
                count += 1
            f.close()
            print("-->There are %d missed jobs" % count)
            if count < 100:
                subprocess.call("bash rerun", shell=True)
                print("-->Submitting %d jobs" % count)
            else:
                print("-->Too many jobs, please check and submit by bash rerun")

        if not os.path.exists("logs"):
            os.makedirs("logs")
        subprocess.call("mv array* logs", shell=True)
        subprocess.call("mv direct_calc_nve* logs", shell=True)
        print("->Wait until remaining NVE trajectories run")
        print("->If this step doesn't work the first time")
        print("->Check your input files")
        open(flag, 'a').close()
        sys.exit()
    return

def step_grabflucts():
    flag='.flag_grabflucts'
    flag_exists = os.path.isfile(flag)
    if flag_exists:
        print("Grab Flucts Flag Exists.")
    else:
        print("Grab Flucts Flag Missing.")
        if not os.path.isfile("flucts.inp"):
            print("flucts.inp missing")
            print("Please generate.")

        subprocess.call("cp %s/src/shell/grabfluctsub.sh ./" % homepath, shell=True)
        os.system("sbatch grabfluctsub.sh")
        print("->Grab flucts is running as a job")
        print("->Please be patient - and check that grab_flucts.py is completed before continuing")
        print("->Grab Flucts has Finished!")
        open(flag, 'a').close()
        sys.exit()
    return

def step_segarray():
    """
    This divides the trajectories into more manageble chunks that speed up things
    like block average calculations and puts them into the SEG directory
    """
    flag='.flag_segarray'
    flag_exists = os.path.isfile(flag)
    if flag_exists:
        print("Segment Creation Array Flag Exists.")
    else:
        print("Segment Creation Array Flag Missing.")
        print("->Running Segement creation calc.")
        if not os.path.exists("SEG"):
            os.makedirs("SEG")
        os.system("sbatch init_segments.sh")
        print("->Segment Creation Array Submitted")
        print("->Please wait until all jobs are complete")
        open(flag, 'a').close()
        sys.exit()
    return

def step_combineseg():
    """
    This combines the segments from the previous step and block averages
    etc to come up with final answers
    """
    flag='.flag_combineseg'
    flag_exists = os.path.isfile(flag)
    if flag_exists:
        print("Segment Combination Flag Exists.")
    else:
        print("Segment Combination Flag Missing.")
        print("->Running segment combination")
        os.system("sbatch combine_segments.sh")
        print("->Segment Combination has finished!")
        open(flag, 'a').close()
        sys.exit()
    return


def step_corrfit():
    flag='.flag_compile'
    flag_exists = os.path.isfile(flag)
    if flag_exists:
        print("Exists.")
    else:
        print("Missing.")


    open(flag, 'a').close()
    return

if __name__ == "__main__":
    mach_name, homepath = step_machinename()
    step_instruct()
    from src.python.read_input import user_input
    inputparam = user_input("input_file")
    step_traj()
    step_nve()
    step_nvecomplete()
    step_checkup()
    step_grabflucts()
    step_segarray()
    step_combineseg()

