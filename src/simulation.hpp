/*
 * File:    simulation.hpp
 * Author:  Allen Sanford (ras9841@rit.edu)
 * Description:
 *      Outlines the Simulation class definition implemented in simulation.cpp
 */

#ifndef __SIMULATION__
#define __SIMULATION__

#include "cells.hpp"
#include "params.hpp"
#include <vector>
#include <random>
#include <stdio.h>

class Simulation {
    std::vector<Cell> cells;
    unsigned int base_seed;
    std::normal_distribution<double> norm_dist;
    std::default_random_engine theta_gen;
    std::default_random_engine phi_gen;
    std::default_random_engine rho_gen;
    void write_cell_loc(FILE *, double);
    void find_collisions(Params::ScalingQuants sq);
    void calc_forces(Params::DimensionlessQuants qd);
    void update_locs(Params::ScalingQuants, Params:: DimensionlessQuants);
public:
    const char *filename;
    // methods
    Simulation(const char *filename, Params::ScalingQuants sq, 
            Params::DimensionlessQuants dq);
    void run(Params::ScalingQuants sq, Params::DimensionlessQuants dq);
};

#endif
