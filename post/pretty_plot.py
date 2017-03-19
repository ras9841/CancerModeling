import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def downsample(name, num):
    time = np.loadtxt(name+"_time.txt")
    hmsd = np.loadtxt(name+"_hmsd.txt")
    cmsd = np.loadtxt(name+"_cmsd.txt")
    hdfc = np.loadtxt(name+"_hdfc.txt")
    cdfc = np.loadtxt(name+"_cdfc.txt")
    
    ntime = []
    ncdfc = []
    nhdfc = []
    ncmsd = []
    nhmsd = []
    count = 0
    while count < time.size:
        if count%num == 0:
            ntime.append(time[count])
            ncdfc.append(cdfc[count])
            nhdfc.append(hdfc[count])
            ncmsd.append(cmsd[count])
            nhmsd.append(hmsd[count])
        count += 1
    
    np.savetxt(name+"_adj_time.txt",ntime,delimiter="\n")
    np.savetxt(name+"_adj_hmsd.txt",nhmsd,delimiter="\n")
    np.savetxt(name+"_adj_cmsd.txt",ncmsd,delimiter="\n")
    np.savetxt(name+"_adj_hdfc.txt",nhdfc,delimiter="\n")
    np.savetxt(name+"_adj_cdfc.txt",ncdfc,delimiter="\n")


def av_dfc():
    name = "long"
    hdfc = np.loadtxt(name+"_hdfc.txt")
    cdfc = np.loadtxt(name+"_cdfc.txt")
    print("mh = %.3f +/- %.3f"%(np.mean(hdfc), np.std(hdfc)))
    print("mc = %.3f +/- %.3f"%(np.mean(cdfc), np.std(cdfc)))

def run():
    bdir = os.path.dirname(__file__)
    base = os.path.join(bdir, '../outputs/')

    N = 3
    tests = ["slowH"]
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
        
        plt.figure(name+" MSD");
        plt.plot(time, h_msd, "-", color=blue)
        plt.plot(time, c_msd, ".-", color=red)
        plt.grid()
        plt.legend(["Healthy Population", "Cancer Population"])
        plt.xlabel(r'Time ($\tau$)')
        plt.ylabel(r'MSD ($\sigma^2$)')

        plt.figure(name+" DFC");
        plt.plot(time, h_dfc, "-", color=blue_l)
        plt.plot(time, c_dfc, ".-", color=red)
        plt.grid()
        plt.legend(["Healthy Population", "Cancer Population"])
        plt.xlabel(r'Time ($\tau$)')
        plt.ylabel(r'MSDFC ($\sigma^2$)')
        plt.show()

if __name__ == "__main__":
    run()
