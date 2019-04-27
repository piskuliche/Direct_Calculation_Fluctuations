#!/usr/bin/env python

import numpy as np
import os,sys,subprocess


def step_machinename():
    exists = os.path.isfile('machine.name')
    if exists:
        f = open('machine.name', 'r')
        mach_name = f.readline().strip()
        f.close()
    else:
        mach_name = str(input("What is MACHINENAME?\n"))
        f = open('machine.name', 'w')
        f.write("%s\n" % mach_name)
        print("Please compile msd_rot_calc, flux_side, and visc_calc")
        sys.exit()
    return

def step_instruct():
    flag='.flag_instruct'
    flag_exists = os.path.isfile(flag)
    if flag_exists:
        print("Instructions Flag exists.")
    else:
        print("Instructions Flag doesn't exist")
        print("Please create input_file")
        print("Please create flucts.inp")
        print("Alternatively, you can run the input_generator to create input_file")
        inpresponse = input("Would you like to do so? y/n\n")
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
        if inputparam.progname == "LAMMPS":
            log_exists = os.path.isfile('log.lammps')
            if log_exists:
                print("Trajectory found")
                open(flag, 'a').close()
            else: 
                print("  At this point, you should run the base trajectory")
                print("  Example NVT and NPT input files are located in src/dependencies")
                sys.exit()
        elif inputparam.progname == "CP2K":
            log_exists = os.path.isfile('cp2k.log')
            if log_exists:
                print("Trajectory found")
                if not os.path.exists("RESTART"):
                    os.makedirs("RESTART")
                subprocess.call("mv *.restart RESTART", shell=True)
                print("  What is the identifier for your restart files?")
                identifier=str(input("  i.e. anything before the # in the .restart files"))
                if os.path.isfile(identifier+str(inputparam.start_config)+".restart"):
                    print("Renaming restart files")
                    for i in range(inputparam.start_config,inputparam.end_config,inputparam.sep_config):
                        os.rename('RESTART/'+identifier+str(i)+'.restart','RESTART/restart.'+str(i))
                    open(flag, 'a').close()
                else:
                    print("Improper identifier, try again.")
                    print("If you would like to skip, please type touch .flag_traj")
            else:
                print("  At this point, you should run the base trajectory")
                print("  Example NVT and NPT input files are located in src/dependencies")
                sys.exit()
        else:
            print("Incorrect program type")
            sys.exit()
    return 

def step_nve():
    flag='.flag_compile'
    flag_exists = os.path.isfile(flag)
    if flag_exists:
        print("Exists.")
    else:
        print("Missing.")


    open(flag, 'a').close()
    return

def step_nvecomplete():
    flag='.flag_compile'
    flag_exists = os.path.isfile(flag)
    if flag_exists:
        print("Exists.")
    else:
        print("Missing.")


    open(flag, 'a').close()
    return

def step_checkup():
    flag='.flag_compile'
    flag_exists = os.path.isfile(flag)
    if flag_exists:
        print("Exists.")
    else:
        print("Missing.")


    open(flag, 'a').close()
    return

def step_grabflucts():
    flag='.flag_compile'
    flag_exists = os.path.isfile(flag)
    if flag_exists:
        print("Exists.")
    else:
        print("Missing.")


    open(flag, 'a').close()
    return

def step_segarray():
    flag='.flag_compile'
    flag_exists = os.path.isfile(flag)
    if flag_exists:
        print("Exists.")
    else:
        print("Missing.")


    open(flag, 'a').close()
    return

def step_combineseg():
    flag='.flag_compile'
    flag_exists = os.path.isfile(flag)
    if flag_exists:
        print("Exists.")
    else:
        print("Missing.")


    open(flag, 'a').close()
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
    step_machinename()
    step_instruct()
    from src.python.read_input import input
    inputparam = input("input_file")
    step_traj()
    step

