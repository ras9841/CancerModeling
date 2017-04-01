#ifndef __TABLE__
#define __TABLE__
#include <vector>       
#include <math.h>
#include "cells.hpp"
#include "simulation.hpp"

class Simulation;

class Table {
    // Attributes
    int num_rho;
    double scale;
    const double pi = std::atan(1)*4.0;
    Simulation *sim;

    // Methods
    void add_cell(Cell *c);
    double atan2(double y, double x);
    double delta_theta(int i);
    double delta_phi(int i);
    int num_theta(int i);
    int num_phi(int i);

    public:
        int num_boxes;
        // Methods
        Table(Simulation *s);
        ~Table();
        std::vector<std::vector<Cell*>> boxes;
        void add_cells();
        void clear_table();
};
#endif
