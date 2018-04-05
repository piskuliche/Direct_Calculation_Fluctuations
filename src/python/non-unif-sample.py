import numpy as np
import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description='''Calculates the log freq''', formatter_class=RawTextHelpFormatter)

parser.add_argument('-start', help = "start")
parser.add_argument('-end', help = "end")
args = parser.parse_args()

start=int(args.start)
end=int(args.end)
end+=50.0
step = 0.0
f=open('time.dat','w')
g=open('real_time.dat','w')
while step < end:
    if step < 500:
        step = step + 10.0
    else:
        step = step + 50.0
    f.write("%d\n" % step)
    g.write("%s\n" % str(step/1000.0))

