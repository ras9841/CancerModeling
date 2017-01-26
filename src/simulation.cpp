#include "simulation.hpp"
#include <time.h>
#include <math.h>

Simulation::Simulation(const char *name, Params::ScalingQuants sq,
                        Params::DimensionlessQuants dq)
{
    // Define constants
    const double pi = std::atan(1)*4;
    double R = dq.Rb;
    #ifdef DEBUG
        R = .5*dq.Rb;
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

void Simulation::write_cell_loc(FILE *file, double time, 
                                Params::ScalingQuants sq) {
    for (unsigned int i=0; i<this->cells.size(); i++) {
        fprintf(file, "%.6f\t%d\t", time*sq.u_time, cells[i].id+1);
        fprintf(file, "%s\t%.6f\t", cells[i].get_type(), cells[i].x*sq.u_length);
        fprintf(file, "%f\t%f\n",cells[i].y*sq.u_length, cells[i].z*sq.u_length);
    }
}


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

void Simulation::calc_forces(Params::DimensionlessQuants dq) {
    for (unsigned int i=0; i<this->cells.size(); i++) {
        cells[i].Fx = 0;
        cells[i].Fy = 0;
        cells[i].Fz = 0;
    }
}

void Simulation::update_locs(Params::ScalingQuants sq, Params::DimensionlessQuants dq) {
    double xnew, ynew, znew;
    for (unsigned int i=0; i<this->cells.size(); i++) {
        #ifdef DEBUG
        if (DEBUG == "brownian") {
            cells[i].dXdt = (sq.u_length*sqrt(2.0/sq.D))*norm_dist(rho_gen);
            cells[i].dYdt = (sq.u_length*sqrt(2.0/sq.D))*norm_dist(rho_gen);
            cells[i].dZdt = (sq.u_length*sqrt(2.0/sq.D))*norm_dist(rho_gen);
        }
        else if (DEBUG == "spp") {
            if (cells[i].is_type(H)) {
                cells[i].dXdt = (dq.prop_H*cos(cells[i].theta)*sin(cells[i].phi))/sqrt(3);
                cells[i].dYdt = (dq.prop_H*sin(cells[i].theta)*sin(cells[i].phi))/sqrt(3);
                cells[i].dZdt = (dq.prop_H*cos(cells[i].phi))/sqrt(3);
            }
            else {
                cells[i].dXdt = (sq.u_length/sq.D)*(dq.prop_C*cos(cells[i].theta)*sin(cells[i].phi))/sqrt(3);
                cells[i].dYdt = (sq.u_length/sq.D)*(dq.prop_C*sin(cells[i].theta)*sin(cells[i].phi))/sqrt(3);
                cells[i].dZdt = (sq.u_length/sq.D)*(dq.prop_C*cos(cells[i].phi))/sqrt(3);
            }
        }
        #else
        // Apply interaction forces
        cells[i].dXdt = (sq.u_length/sq.u_energy)*cells[i].Fx;
        cells[i].dYdt = (sq.u_length/sq.u_energy)*cells[i].Fy;
        cells[i].dZdt = (sq.u_length/sq.u_energy)*cells[i].Fz;
        // Apply self propulsion
        if (cells[i].is_type(H)) {
            cells[i].dXdt += (dq.prop_H*cos(cells[i].theta)*sin(cells[i].phi))/sqrt(3);
            cells[i].dYdt += (dq.prop_H*sin(cells[i].theta)*sin(cells[i].phi))/sqrt(3);
            cells[i].dZdt += (dq.prop_H*cos(cells[i].phi))/sqrt(3);
        }
        else {
            cells[i].dXdt += (sq.u_length/sq.D)*(sq.u_length/sq.D)*(dq.prop_C*cos(cells[i].theta)*sin(cells[i].phi))/sqrt(3);
            cells[i].dYdt += (sq.u_length/sq.D)*(dq.prop_C*sin(cells[i].theta)*sin(cells[i].phi))/sqrt(3);
            cells[i].dZdt += (sq.u_length/sq.D)*(dq.prop_C*cos(cells[i].phi))/sqrt(3);
        }
        // Apply diffusion
        cells[i].dXdt += (sq.u_length*sqrt(2.0/sq.D))*norm_dist(rho_gen);
        cells[i].dYdt += (sq.u_length*sqrt(2.0/sq.D))*norm_dist(rho_gen);
        cells[i].dZdt += (sq.u_length*sqrt(2.0/sq.D))*norm_dist(rho_gen);
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
        cells[i].theta += dq.dt*(sq.u_length*sqrt(6/sq.D))*norm_dist(theta_gen);
        cells[i].phi += dq.dt*(sq.u_length*sqrt(6/sq.D))*norm_dist(theta_gen);
    }
}


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
        write_cell_loc(file, t, sq);        
        find_collisions(sq);
        calc_forces(dq);
        update_locs(sq, dq);
        t += dq.dt;
    }
    write_cell_loc(file, t, sq);        

    fclose(file);
}



