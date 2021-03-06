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
#include "table.hpp"
#include <vector>
#include <random>
#include <stdio.h>
#include <math.h>

class Table; 

class Simulation {
    // attributes
    unsigned int base_seed;
    std::normal_distribution<double> norm_dist;
    std::default_random_engine theta_gen;
    std::default_random_engine phi_gen;
    std::default_random_engine rho_gen;
    const double pi = std::atan(1)*4.0;
    int pop_size;

    // methods
    void write_cell_loc(FILE *, double);
    void write_cell_msd(FILE *, double);
    void write_cell_dfc(FILE *, double);
    void find_collisions(Table *t);
    void calc_forces();
    void update_locs();
    void divide_cells(CellType T);
public:
    // attributes
    const char *loc_name;
    const char *msd_name;
    const char *dfc_name;
    Params::SysQuants sq; 
    Params::DimensionlessQuants dq;
    std::vector<Cell> cells;

    // methods
    Simulation(){};
    Simulation(const char *name_loc, const char *name_msd, 
            const char *name_dfc, Params::SysQuants sq, 
            Params::DimensionlessQuants dq);
    void run();
};

#endif
