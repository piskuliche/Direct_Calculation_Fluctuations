import numpy as np
import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description='''Calculates the log freq''', formatter_class=RawTextHelpFormatter)

parser.add_argument('-start', help = "start")
parser.add_argument('-end', help = "end")
parser.add_argument('-step', help = "step")
args = parser.parse_args()

start=int(args.start)
end=int(args.end)
dt=float(args.step)

step = 0.0
f=open('time.dat','w')
g=open('real_time.dat','w')
while step < end:
    f.write("%d\n" % step)
    g.write("%s %s\n" % (str(step/1000.0),1.0))
    step = step + dt*1000.0
