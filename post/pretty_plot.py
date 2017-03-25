import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def compute_velocity(time, dfc):
    disp = np.sqrt(dfc)
    vel = np.zeros(time.size)
    for i in range(1,time.size):
        vel[i] = (disp[i]-disp[i-1])/(time[i]-time[i-1])
    return vel

def run():
    bdir = os.path.dirname(__file__)
    base = os.path.join(bdir, '../outputs/')

    N = 1
    tests = ["many"]
    #tests = ["prop", "elast", "surf", "slow", "long_adj"]
    
    for name in tests:
        msd = np.loadtxt(base+name+"-%d_msd.csv"%(1), comments="#")
        time = msd[:,0]
        h_msd = np.zeros([time.size, 1])
        c_msd = np.zeros([time.size, 1])
        h_dfc = np.zeros([time.size, 1])
        c_dfc = np.zeros([time.size, 1])
        
        for i in range(1,N+1):
            msd = np.loadtxt(base+name+"-%d_msd.csv"%(i), comments="#")
            h_msd = h_msd + msd[:, 1]
            c_msd = c_msd + msd[:, 2]
            dfc = np.loadtxt(base+name+"-%d_dfc.csv"%(i), comments="#")
            h_dfc = h_dfc + dfc[:, 1]
            c_dfc = c_dfc + dfc[:, 2]
        
        h_msd = np.transpose(h_msd/N)
        c_msd = np.transpose(c_msd/N)
        h_dfc = np.transpose(h_dfc/N)
        c_dfc = np.transpose(c_dfc/N)
        
        mpl.rcParams["font.size"] = 18
        blue = (0, 0, 153./256, .75)
        red = (204./256, 0, 0, .75)
        blue_l = (0, 0, 153./256, .3)
        
        time = np.array([t for t in time])
        h_msd = [h[0] for h in h_msd]
        c_msd = [c[0] for c in c_msd]
        h_dfc = [h[0] for h in h_dfc]
        c_dfc = [c[0] for c in c_dfc]
        
        plt.figure(name+" MSD");
        plt.plot(time, h_msd, "-", color=blue)
        plt.plot(time, c_msd, ".-", color=red)
        plt.grid()
        plt.legend(["Healthy Population", "Cancer Population"],\
               numpoints=1) 
        plt.xlabel(r'Time ($\tau$)')
        plt.ylabel(r'MSD ($\sigma^2$)')

        plt.figure(name+" DFC");
        plt.plot(time, h_dfc, "-", color=blue_l)
        plt.plot(time, c_dfc, ".-", color=red)
        plt.grid()
        plt.legend(["Healthy Population", "Cancer Population"])
        plt.xlabel(r'Time ($\tau$)')
        plt.ylabel(r'MSDFC ($\sigma^2$)')

        plt.figure(name+" Velocity");
        plt.plot(time, compute_velocity(time,h_dfc), "-", color=blue_l)
        plt.plot(time, compute_velocity(time,c_dfc), ".-", color=red)
        plt.grid()
        plt.legend(["Healthy Population", "Cancer Population"])
        plt.xlabel(r'Time ($\tau$)')
        plt.ylabel(r'Velocity ($\sigma/\tau$)')
        plt.show()
if __name__ == "__main__":
    run()
