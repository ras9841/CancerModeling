import plot_trajectories as traj
import pretty_plot as plots

if __name__ == "__main__":
    tests = ["elastic_division"]
    #tests = ["diff_prop", "diff_elast", "diff_surf_E", "full"]
    num_tests = 1
    size_per_pop = 128
    pac_frac = 0.9

    plots.run(num_tests, tests)
    #traj.run(size_per_pop, pac_frac, tests)
