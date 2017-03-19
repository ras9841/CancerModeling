import os
import matplotlib.pyplot as plt

def read_data(test):
    hlocs = [[] for _ in range(128)]
    clocs = [[] for _ in range(128)]
    time = {}

    bdir = os.path.dirname(__file__)
    base = os.path.join(bdir, '../outputs/')
    with open(base+test+"_loc.csv","r") as f:
        for _ in range(4): # skip header
            f.readline()
        for line in f:
            t, cid, ctype, x, y, z = line.strip().split("\t")
            if t not in time.keys():
                time[t] = 1
            if ctype == "H":
                hlocs[int(cid)-1].append((float(x), float(y), float(z)))
            else:
                clocs[int(cid)-1-128].append((float(x), float(y), float(z)))
        return time.keys, hlocs, clocs
if __name__ == "__main__":
    time, hlocs, clocs = read_data("full-1")
    cell_num = 128
