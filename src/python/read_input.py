class input:

    def __init__(self,inputfile):
        ip = open(inputfile)

        self.molec=[]
        for i, line in enumerate(ip):
            if i == 6:
                self.prog=str(line.strip())
            if i == 10:
                # 11th line
                self.num_molecs=int(line)
            if i == 12:
                self.molec.append(str(line.strip()))
            if i == 13 and self.num_molecs > 1:
                self.molec.append(str(line.strip()))
            if i == 14 and self.num_molecs > 2:
                self.molec.append(str(line.strip()))
            if i == 15 and self.num_molecs > 3:
                self.molec.append(str(line.strip()))
            if i == 17:
                self.start_config=int(line.strip())
            if i == 19:
                self.end_config=int(line.strip())
            if i == 21:
                self.sep_config=float(line.strip())
            if i == 23:
                self.timestep=float(line.strip())
            if i == 25:
                self.num_times=int(line.strip())
            if i == 27:
                self.nblocks=int(line.strip())
            if i == 29:
                self.num_files=int(line.strip())
            if i == 31:
                self.nve_length=int(line.strip())
        ip.close()
