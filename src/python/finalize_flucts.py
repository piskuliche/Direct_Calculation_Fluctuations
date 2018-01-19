import numpy as np


f = open("molnames.txt","r") #opens file with name of "test.txt"

molnames = []

for line in f:
        print "python flucts_calc.py -inp test.inp -files 5000 -blocks 10 -mol %s\n" % line.strip()
