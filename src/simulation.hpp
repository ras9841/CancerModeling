#ifndef __SIMULATION__
#define __SIMULATION__

#include "cells.hpp"
#include "params.hpp"
#include <vector>
#include <stdio.h>

class Simulation {
    std::vector<Cell> cells;
    unsigned int base_seed;
    void write_cell_loc(FILE *, double);
    void find_collisions(Params::ScalingQuants sq);
    void calc_forces();
    void update_locs();
public:
    const char *filename;
    // methods
    Simulation(const char *filename, Params::ScalingQuants sq, 
            Params::DimensionlessQuants dq);
    void run(Params::ScalingQuants sq, Params::DimensionlessQuants dq);
};

#endif
