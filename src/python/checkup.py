import numpy as np
import os.path

with open("file_names") as f:
    for file in f:
        filepath=str("FILES/")+str(file).strip()+str("/traj_")+str(file).strip()+str(".xyz")
        if os.path.isfile(filepath):
            output=1
        else:
            cd=str("cd FILES/")+str(file).strip()
            output=str("msub ")+str("water_nve.sh")
            back="cd ../../"
            print cd
            print output
            print back
        #DONE
    #DONE

#DONE
