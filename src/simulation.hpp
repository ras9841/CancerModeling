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
#include <math.h>

class Simulation {
    // attributes
    std::vector<Cell> cells;
    unsigned int base_seed;
    std::normal_distribution<double> norm_dist;
    std::default_random_engine theta_gen;
    std::default_random_engine phi_gen;
    std::default_random_engine rho_gen;
    Params::SysQuants sq; 
    Params::DimensionlessQuants dq;
    const double pi = std::atan(1)*4.0;

    // methods
    void write_cell_loc(FILE *, double);
    void find_collisions();
    void calc_forces();
    void update_locs();
public:
    // attributes
    const char *filename;
    // methods
    Simulation(const char *filename, Params::SysQuants sq, 
            Params::DimensionlessQuants dq);
    void run();
};

#endif
