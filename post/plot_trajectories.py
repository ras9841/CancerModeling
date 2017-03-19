import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np 
import math

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
                time[t] = float(t)
            if ctype == "H":
                hlocs[int(cid)-1].append((float(x), float(y), float(z)))
            else:
                clocs[int(cid)-1-128].append((float(x), float(y), float(z)))
        return time.values(), hlocs, clocs

def plot_cell_traj(ax, time, num, cells, sty):
    xvals = []
    yvals = []
    zvals = []

    for i in range(len(time)):
        xvals.append(cells[num][i][0])
        yvals.append(cells[num][i][1])
        zvals.append(cells[num][i][2])

    ax.scatter(xvals[0], yvals[0], zvals[0], color="g", s=100)
    ax.plot(xvals, yvals, zvals, sty)
    ax.scatter(xvals[-1], yvals[len(time)-1], zvals[-1], color="r", s=100)
    return ax

def plot_sphere(ax, R):
    x = []
    y = []
    z = []
    for theta in np.linspace(0, 2*math.pi, 50):
        for phi in np.linspace(0, math.pi, 50):
            x.append(R*math.cos(theta)*math.sin(phi))
            y.append(R*math.sin(theta)*math.sin(phi))
            z.append(R*math.cos(phi))
    ax.plot_wireframe(x,y,z, color="k", alpha=.3)

def run():
    R  = 8.3
    time, hlocs, clocs = read_data("slowH-1")
    cell_num = 0
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.set_aspect("equal")
    ax = plot_cell_traj(ax, list(time), 0, clocs, "m-")
    ax = plot_cell_traj(ax, list(time), 5, clocs, "m-")
    ax = plot_cell_traj(ax, list(time), 10, clocs, "m-")
    #ax = plot_cell_traj(ax, list(time), 0, hlocs, "b-")
    #ax = plot_cell_traj(ax, list(time), 5, hlocs, "b-")
    #ax = plot_cell_traj(ax, list(time), 10, hlocs, "b-")
    ax = plot_sphere(ax, R)
    plt.show()

if __name__ == "__main__":
    run()
