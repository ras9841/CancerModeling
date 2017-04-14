import plot_trajectories as traj
import pretty_plot as plots

if __name__ == "__main__":
    #tests = ["full"] 
    tests = ["division"] 
    #tests = ["small_HC", "big_HC"] 
    #tests = ["diff_surf_E", "full"] 
    #tests = ["diff_prop", "diff_elast", "diff_surf_E", "full"]
    num_tests = 2
    size_per_pop = 100
    pac_frac = 0.95

    plots.run(num_tests, tests)
    traj.run(size_per_pop, pac_frac, tests)
