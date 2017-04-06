import plot_trajectories as traj
import pretty_plot as plots

if __name__ == "__main__":
    tests = ["full", "division"] 
    #tests = ["diff_prop", "diff_elast", "diff_surf_E", "full"]
    num_tests = 10
    size_per_pop = 1000
    pac_frac = 0.95

    plots.run(num_tests, tests)
    #traj.run(size_per_pop, pac_frac, tests)
