#include "simulation.hpp"
#include <time.h>
#include <math.h>
#include <random>

Simulation::Simulation(const char *name, Params::ScalingQuants sq,
                        Params::DimensionlessQuants dq)
{
    // Define pi
    const double pi = std::atan(1)*4;

    // Initialize attributes
    this->filename = name;

    // Setup Simulation
    unsigned int base = time(NULL);
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dist (0, 1);

    this->base_seed = base;

    // Generate rho values
    double rho_vals[2*sq.N];
    
    gen.seed(base); // H
    for (int i = 0; i<sq.N; i++) { 
        rho_vals[i] = dq.Rb*cbrt(dist(gen)); 
    } 
    gen.seed(base+1); // C
    for (int i = sq.N; i<2*sq.N; i++) { 
        rho_vals[i] = dq.Rb*cbrt(dist(gen)); 
    } 
    
    // Generate theta values
    double theta_vals[2*sq.N];
    
    gen.seed(base+2); // H
    for (int i = 0; i<sq.N; i++) { 
        theta_vals[i] = 2*pi*dist(gen); 
    } 
    gen.seed(base+3); // C
    for (int i = sq.N; i<2*sq.N; i++) { 
        theta_vals[i] = 2*pi*dist(gen); 
    } 
    
    // Generate phi values
    double phi_vals[2*sq.N];
    
    gen.seed(base+4); // H
    for (int i = 0; i<sq.N; i++) { 
        phi_vals[i] = acos(1-2*dist(gen)); 
    } 
    gen.seed(base+5); // C
    for (int i = sq.N; i<2*sq.N; i++) { 
        phi_vals[i] = acos(1-2*dist(gen)); 
    } 

    double x, y, z;
    for (int i = 0; i<sq.N; i++) {
        x = rho_vals[i] * cos(theta_vals[i]) * sin(phi_vals[i]);
        y = rho_vals[i] * sin(theta_vals[i]) * sin(phi_vals[i]);
        z = rho_vals[i] * cos(phi_vals[i]);
        this->cells.push_back(Cell(i, H, x, y, z));
    }
    for (int i = sq.N; i<2*sq.N; i++) {
        x = rho_vals[i] * cos(theta_vals[i]) * sin(phi_vals[i]);
        y = rho_vals[i] * sin(theta_vals[i]) * sin(phi_vals[i]);
        z = rho_vals[i] * cos(phi_vals[i]);
        this->cells.push_back(Cell(i+1, C, x, y, z));
    }
}

void Simulation::write_cell_loc(FILE *file, double time) {
    for (unsigned int i=0; i<this->cells.size(); i++) {
        fprintf(file, "%.6f\t%d\t%s\t", time, cells[i].id+1, cells[i].get_type());
        fprintf(file, "%.6f\t%f\t%f\n", cells[i].x, cells[i].y, cells[i].z);
    }
}

/*
void find_collisions(Params::ScalingQuants dq) {
    for (unsigned int i=0; i<this->cells.size(); i++) {
        cell[i].adjlst.clear();
        for (unsigned int j=0; j<this->cells.size(); j++) {
            if (i != j) {
                if (cell[i].hits(cell[j])) {
                    cell[i].addjlst.push_back(cell[j].id);
                }
            }
    }
}

void calc_forces();
void update_locs();
*/

void Simulation::run(Params::ScalingQuants sq, Params::DimensionlessQuants dq) {
    // Get current time
    time_t rawtime;
    struct tm * timeinfo;
    char buffer [80];

    time (&rawtime);
    timeinfo = localtime (&rawtime);
    strftime (buffer,80,"%c",timeinfo);
    
    FILE *file = fopen(this->filename, "w");
    fprintf(file, "# Time started: %s\n", buffer); 
    fprintf(file, "# Simulation using %d cells\n", 2*sq.N);
    fprintf(file, "# Format: time(tau)\t cellID\t cellType\t x\t y\t z\n");

    double t = 0;
    while (t < 10*dq.dt) {
        write_cell_loc(file, t);        
        //find_collisions();
        //calc_forces();
        //update_locs();
        t += dq.dt;
    }
    //write_cell_loc(file, t);        

    fclose(file);
}



