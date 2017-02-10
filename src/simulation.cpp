/*
 * File:    simulation.cpp
 * Author:  Allen Sanford (ras9841@rit.edu)
 * Description:
 *      Defines the simulation class and its methods.
 *      Contains all logic for the simlutaion.
 */

// Imports
#include "simulation.hpp"
#include <time.h>
#include <math.h>

/// Constructor for a Simulation objection.
///
/// The initial setup for the simulation is created using the values specified 
/// by the scaling and dimensionless quantities. The cells are scattered 
/// uniformly throughout the volume of bounding sphere. The coordinates are 
/// initially generated using spherical coordinates (each with an independent 
/// distribution) but are saved in the cells in cartesian coordinates in order 
/// to match the force model. All random generation is based off a base seed, 
/// which is generated based on the computing system's time.
///
/// Keyword Arguments:
///     name    --  name used for output files.
///     sq      --  scaling quantities used for length, time, and energy scales.
///     dq      --  dimensionless quantities used in the simulation
Simulation::Simulation(const char *name, Params::ScalingQuants sq,
                        Params::DimensionlessQuants dq)
{
    // Define constants
    const double pi = std::atan(1)*4;
    double R = dq.Rb;
    #ifdef DEBUG
        R = .2*dq.Rb;
    #endif

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
        rho_vals[i] = R*cbrt(dist(gen)); 
    } 
    gen.seed(base+1); // C
    for (int i = sq.N; i<2*sq.N; i++) { 
        rho_vals[i] = R*cbrt(dist(gen)); 
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

    // Setup distributions for vp orientations
    this->theta_gen.seed(base+6);
    this->phi_gen.seed(base+7);
    this->rho_gen.seed(base+8);

    double x, y, z;
    for (int i = 0; i<sq.N; i++) {
        x = rho_vals[i] * cos(theta_vals[i]) * sin(phi_vals[i]);
        y = rho_vals[i] * sin(theta_vals[i]) * sin(phi_vals[i]);
        z = rho_vals[i] * cos(phi_vals[i]);
        this->cells.push_back(Cell(i, H, x, y, z, 
                    norm_dist(theta_gen), norm_dist(phi_gen)));
    }
    for (int i = sq.N; i<2*sq.N; i++) {
        x = rho_vals[i] * cos(theta_vals[i]) * sin(phi_vals[i]);
        y = rho_vals[i] * sin(theta_vals[i]) * sin(phi_vals[i]);
        z = rho_vals[i] * cos(phi_vals[i]);
        this->cells.push_back(Cell(i+1, C, x, y, z,
                    norm_dist(theta_gen), norm_dist(phi_gen)));
    }
}

/// Updates the results file with the current state of the system.
///
/// The output quantities (time,x,y,z) are with respect to the scaling
/// quantities ulength and utime.
///
/// Keyword Arguments:
///     file    --  file pointer to the output data location.
///     time    --  current time in units of tau
void Simulation::write_cell_loc(FILE *file, double time) { 
    for (unsigned int i=0; i<this->cells.size(); i++) {
        fprintf(file, "%.6f\t%d\t", time, cells[i].id+1);
        fprintf(file, "%s\t%.6f\t", cells[i].get_type(), cells[i].x);
        fprintf(file, "%f\t%f\n", cells[i].y, cells[i].z);
    }
}

/// Finds all cells in the current simulation state that are colliding.
///
/// Current implementation in a brute force O(N^2) direct comparison. Collided 
/// cells are added to the cells' adjacency list for future use in the force
/// calculations.
///
/// Keyword Arguments:
///     sq  --  scaling quantities (used for unit length)
void Simulation::find_collisions(Params::ScalingQuants sq) {
    for (unsigned int i=0; i<this->cells.size(); i++) {
        cells[i].adjlst.clear();
        for (unsigned int j=0; j<this->cells.size(); j++) {
            if (i != j) {
                if (cells[i].hits(cells[j], sq.u_length)) {
                    cells[i].adjlst.push_back(cells[j].id);
                }
            }
        }
    }
}

/// Calculates the adhesive and repulsive forces for each cell.
///
/// Currently, these interactions are set to zero. Future versions will include
/// the JKR model for adhesion and repulsion.
///
/// Keywork Arguments:
///     dq  --  demensionless quantities used in force calculation.
void Simulation::calc_forces(Params::DimensionlessQuants dq) {
    for (unsigned int i=0; i<this->cells.size(); i++) {
        cells[i].Fx = 0;
        cells[i].Fy = 0;
        cells[i].Fz = 0;
    }
}

/// Update cell positions by integrating the 3D Langevin equations.
///
/// Currently, forces contributing to the particles velocity include Brownian 
/// motion, self-propulsion, cell-cell adhession, and cell-cell repulsion. The
/// numerical integration is done with Forward Euler O(dt).
///
/// Keyword Arguments:
///     sq  --  scaling quantity used to nondimensionalize the force equations
///     dq  --  demensionless quantities used in force calculation.
void Simulation::update_locs(Params::ScalingQuants sq, Params::DimensionlessQuants dq) {
    double xnew, ynew, znew;
    for (unsigned int i=0; i<this->cells.size(); i++) {
        #ifdef DEBUG
        if (DEBUG == "brownian") {
            cells[i].dXdt = (sq.u_length*sqrt(2.0/(sq.D*dq.dt*sq.u_time)))*norm_dist(rho_gen);
            cells[i].dYdt = (sq.u_length*sqrt(2.0/(sq.D*dq.dt*sq.u_time)))*norm_dist(rho_gen);
            cells[i].dZdt = (sq.u_length*sqrt(2.0/(sq.D*dq.dt*sq.u_time)))*norm_dist(rho_gen);
        }
        else if (DEBUG == "spp") {
            if (cells[i].is_type(H)) {
                cells[i].dXdt = (sq.u_length/sq.D)*(dq.prop_H*cos(cells[i].theta)*sin(cells[i].phi))/sqrt(3);
                cells[i].dYdt = (sq.u_length/sq.D)*(dq.prop_H*sin(cells[i].theta)*sin(cells[i].phi))/sqrt(3);
                cells[i].dZdt = (sq.u_length/sq.D)*(dq.prop_H*cos(cells[i].phi))/sqrt(3);
            }
            else {
                cells[i].dXdt = (sq.u_length/sq.D)*(dq.prop_C*cos(cells[i].theta)*sin(cells[i].phi))/sqrt(3);
                cells[i].dYdt = (sq.u_length/sq.D)*(dq.prop_C*sin(cells[i].theta)*sin(cells[i].phi))/sqrt(3);
                cells[i].dZdt = (sq.u_length/sq.D)*(dq.prop_C*cos(cells[i].phi))/sqrt(3);
            }
        }
        #else
        cells[i].dXdt = (sq.u_length/sq.u_energy)*cells[i].Fx;
        cells[i].dYdt = (sq.u_length/sq.u_energy)*cells[i].Fy;
        cells[i].dZdt = (sq.u_length/sq.u_energy)*cells[i].Fz;
        // Apply self propulsion
        if (cells[i].is_type(H)) {
            cells[i].dXdt = (sq.u_length/sq.D)*(dq.prop_H*cos(cells[i].theta)*sin(cells[i].phi))/sqrt(3);
            cells[i].dYdt = (sq.u_length/sq.D)*(dq.prop_H*sin(cells[i].theta)*sin(cells[i].phi))/sqrt(3);
            cells[i].dZdt = (sq.u_length/sq.D)*(dq.prop_H*cos(cells[i].phi))/sqrt(3);
        }
        else {
            cells[i].dXdt = (sq.u_length/sq.D)*(dq.prop_C*cos(cells[i].theta)*sin(cells[i].phi))/sqrt(3);
            cells[i].dYdt = (sq.u_length/sq.D)*(dq.prop_C*sin(cells[i].theta)*sin(cells[i].phi))/sqrt(3);
            cells[i].dZdt = (sq.u_length/sq.D)*(dq.prop_C*cos(cells[i].phi))/sqrt(3);
        }
        // Apply diffusion
        cells[i].dXdt = (sq.u_length*sqrt(2.0/(sq.D*dq.dt*sq.u_time)))*norm_dist(rho_gen);
        cells[i].dYdt = (sq.u_length*sqrt(2.0/(sq.D*dq.dt*sq.u_time)))*norm_dist(rho_gen);
        cells[i].dZdt = (sq.u_length*sqrt(2.0/(sq.D*dq.dt*sq.u_time)))*norm_dist(rho_gen);
        #endif
        
        // Update positions and orientations (forward euler)
        xnew = cells[i].x + dq.dt*cells[i].dXdt;
        ynew = cells[i].y + dq.dt*cells[i].dYdt;
        znew = cells[i].z + dq.dt*cells[i].dZdt;

        // Validity check: stay in sphere
        if (pow(xnew,2)+pow(ynew,2)+pow(znew,2) < pow(dq.Rb,2))
        {
            cells[i].x = cells[i].x + dq.dt*cells[i].dXdt;
            cells[i].y = cells[i].y + dq.dt*cells[i].dYdt;
            cells[i].z = cells[i].z + dq.dt*cells[i].dZdt;
        }

        // Update theta and phi
        cells[i].theta += dq.dt*(sq.u_length*sqrt(6/(sq.D*dq.dt*sq.u_time)))*norm_dist(theta_gen);
        cells[i].phi += dq.dt*(sq.u_length*sqrt(6/(sq.D*dq.dt*sq.u_time)))*norm_dist(theta_gen);    
    }
}

/// Logic for running the simulation
/// 
/// Creates the filename used to write out particle position. This includes 
/// opening and closing the file pointer. Order of operations: (1) write out
/// current cell locations, (2) find cells that are colliding, (3) calculate
/// the adhesive and repulsive forces acting on each cell, and (4) update
/// the particle locations.
///
/// Keyword Arguments:
///     sq  --  scaling quantity used to nondimensionalize the force equations
///     dq  --  demensionless quantities used in force calculation.
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
    fprintf(file, "# Bounding Radius: %f\t Time Scale: %f\n", dq.Rb, sq.u_time);
    fprintf(file, "# Format: time(tau)\t cellID\t cellType\t x\t y\t z\n");

    double t = 0;
    while (t < dq.tf) {
        write_cell_loc(file, t);        
        find_collisions(sq);
        calc_forces(dq);
        update_locs(sq, dq);
        t += dq.dt;
    }
    write_cell_loc(file, t);        

    fclose(file);
}



