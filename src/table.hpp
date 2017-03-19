#ifndef __TABLE__
#define __TABLE__
#include "simulation.hpp"
#include "linked_list.hpp"
#include <math.h>

class Simulation::Table {
    // Attributes
    int num_rho;
    double scale;
    const double pi = std::atan(1)*4.0;
    
    // Methods
    void add_cell(Cell *c);
    double atan2(double y, double x);
    double delta_theta(int i);
    double delta_phi(int i);
    int num_theta(int i);
    int num_phi(int i);

    public:
        int num_boxes;
        Simulation *sim;
        // Methods
        Table(Simulation *s);
        ~Table();
        LinkedList<Cell*> *boxes;
        void add_cells();
        void clear_table();
};
#endif
