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

def pretty_plot():
    #tests = ["long_adj"]
    tests = ["prop", "elast", "surf", "slow", "long_adj"]
    for name in tests:
        time = np.loadtxt(name+"_time.txt")
        hmsd = np.loadtxt(name+"_hmsd.txt")
        cmsd = np.loadtxt(name+"_cmsd.txt")
        hdfc = np.loadtxt(name+"_hdfc.txt")
        cdfc = np.loadtxt(name+"_cdfc.txt")

        mpl.rcParams["font.size"] = 18
        blue = (0, 0, 153./256, .75)
        red = (204./256, 0, 0, .75)
        blue_l = (0, 0, 153./256, .3)
        
        plt.figure(name+" MSD");
        plt.plot(time, hmsd, "-", color=blue,label="Healthy Population")
        plt.plot(time, cmsd, ".-", color=red, label="Cancer Population")
        plt.grid()
        plt.legend()
        plt.xlabel(r'Time ($\tau$)')
        plt.ylabel(r'<||r(t)-r(0)||> ($\sigma^2$)')

        plt.figure(name+" DFC");
        plt.plot(time, hdfc, "-", color=blue_l,label="Healthy Population")
        plt.plot(time, cdfc, ".-", color=red, label="Cancer Population")
        plt.grid()
        plt.legend()
        plt.xlabel(r'Time ($\tau$)')
        plt.ylabel(r'<||r(t)||> ($\sigma^2$)')
        plt.show()


def plot_rad():
    name = "slow2"
    time = np.loadtxt(name+"_time.txt")
    print(time)

if __name__ == "__main__":
    #downsample("long", 3)
    #pretty_plot()
    #av_dfc()
    plot_rad()
