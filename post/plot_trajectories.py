import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np 
import math

def read_data(N, test):
    hlocs = [[] for _ in range(3)]
    clocs = [[] for _ in range(3)]

    bdir = os.path.dirname(__file__)
    base = os.path.join(bdir, '../outputs/')
    time = {} 
    
    with open(base+test+"_loc.csv","r") as f:
        for _ in range(4): # skip header
            f.readline()
        for line in f:
            t, cid, ctype, x, y, z = line.strip().split("\t")
            if t not in time.keys():
                time[t] = float(t)
    times = list(time.values())
    times.sort()
    times = [times[0], times[int(len(times)/2)] , times[-1]]
    time_dict = {}
    i = 0
    for t in times:
        time_dict[t] = i
        i += 1

    with open(base+test+"_loc.csv","r") as f:
        for _ in range(4): # skip header
            f.readline()
        for line in f:
            t, cid, ctype, x, y, z = line.strip().split("\t")
            if float(t) in times:
                i = time_dict[float(t)] 
                if ctype == "H":
                    hlocs[i].append([float(x), float(y), float(z)])
                else:
                    clocs[i].append([float(x), float(y), float(z)])
    return times, hlocs, clocs

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
    return ax

def plot_cells(ax, cells, color):
    get_x = lambda x: x[0]
    get_y = lambda x: x[1]
    get_z = lambda x: x[2]
    xvals = list(map(get_x, cells))
    yvals = list(map(get_y, cells))
    zvals = list(map(get_z, cells))
    ax.scatter(xvals, yvals, zvals, color=color, s=50)
    return ax

def run(N, pac_frac, test_names): 
    blue = (0,0,153./250., .80)
    red = (204./256.,0,0., .80)

    for name in test_names:
        times, hlocs, clocs = read_data(N, name+"-1")
        R = (4*N/pac_frac)**(1./3.)+.5

        for i in range(3):
            fig = plt.figure()
            ax = fig.add_subplot(111, projection="3d")
            ax.set_aspect("equal")
            ax = plot_sphere(ax, R)
            ax = plot_cells(ax, hlocs[i], blue)
            ax = plot_cells(ax, clocs[i], red)
            ax.set_title(r'Locations at t=%d$\tau$'%(times[i]))
            #plt.legend(['Bounding Sphere', 'Healthy Cell', 'Cancer Cell'],
            #        loc=5)
            plt.grid('off')
            plt.axis('off')
            plt.show()

if __name__ == "__main__":
    #run(128, .9, ["slow_hash"])
    run(1000, .9, ["many"])
