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
///     name_loc    --  name used for cell location output file.
///     name_msd    --  name used for cell MSD output file.
///     name_dfc    --  name used for cell DFC output file.
///     sq      --  system quantities used for length, time, and energy scales.
///     dq      --  dimensionless quantities used in the simulation
Simulation::Simulation(const char *name_loc, const char *name_msd, 
                    const char *name_dfc, Params::SysQuants sq,
                    Params::DimensionlessQuants dq)
{
    // Define constants
    double R = dq.Rb;
    #ifdef DEBUG
        if (DEBUG == "brownian" || DEBUG == "spp") {R = .2*dq.Rb;}
    #endif

    // Initialize attributes
    this->loc_name = name_loc;
    this->msd_name = name_msd;
    this->dfc_name = name_dfc;
    this->sq = sq;
    this->dq = dq;
    this->pop_size = 2*sq.N;

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
    //gen.seed(base+1); // C
    for (int i = sq.N; i<2*sq.N; i++) { 
        rho_vals[i] = R*cbrt(dist(gen)); 
    } 
    
    // Generate theta values
    double theta_vals[2*sq.N];
    
    //gen.seed(base+2); // H
    for (int i = 0; i<sq.N; i++) { 
        theta_vals[i] = 2*pi*dist(gen); 
    } 
    //gen.seed(base+3); // C
    for (int i = sq.N; i<2*sq.N; i++) { 
        theta_vals[i] = 2*pi*dist(gen); 
    } 
    
    // Generate phi values
    double phi_vals[2*sq.N];
    
    //gen.seed(base+4); // H
    for (int i = 0; i<sq.N; i++) { 
        phi_vals[i] = acos(1-2*dist(gen)); 
    } 
    //gen.seed(base+5); // C
    for (int i = sq.N; i<2*sq.N; i++) { 
        phi_vals[i] = acos(1-2*dist(gen)); 
    } 

    // Setup distributions for vp orientations
    this->theta_gen.seed(base+1);
    this->phi_gen.seed(base+2);
    this->rho_gen.seed(base+3);

    double x, y, z;
    for (int i = 0; i<sq.N; i++) { // Healthy
        x = rho_vals[i] * cos(theta_vals[i]) * sin(phi_vals[i]);
        y = rho_vals[i] * sin(theta_vals[i]) * sin(phi_vals[i]);
        z = rho_vals[i] * cos(phi_vals[i]);
        this->cells.push_back(Cell(i, H, x, y, z, 
                    norm_dist(theta_gen), norm_dist(phi_gen), sq.radius_H,
                    sq.e_mod_H, sq.poisson_H, sq.prop_H 
                    ));
    }
    for (int i = sq.N; i<2*sq.N; i++) { // Cancer
        x = rho_vals[i] * cos(theta_vals[i]) * sin(phi_vals[i]);
        y = rho_vals[i] * sin(theta_vals[i]) * sin(phi_vals[i]);
        z = rho_vals[i] * cos(phi_vals[i]);
        this->cells.push_back(Cell(i, C, x, y, z,
                    norm_dist(theta_gen), norm_dist(phi_gen), sq.radius_C,
                    sq.e_mod_C, sq.poisson_C, sq.prop_C 
                    ));
    }
}

/// Updates the results file with the current state of the system.
///
/// The output quantities (time,id,type,x,y,z) with respect to the scaling
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
/// Updates the results file with the current state of the system.
///
/// The output quantities (time, msd_H, msd_C) with respect to the scaling
/// quantities ulength and utime.
///
/// Keyword Arguments:
///     file    --  file pointer to the output data location.
///     time    --  current time in units of tau
void Simulation::write_cell_msd(FILE *file, double time) { 
    double count_H = 0.0, count_C = 0.0;
    double sum_msd_H = 0.0, sum_msd_C = 0.0;
    double msd = 0.0;
    for (unsigned int i=0; i<this->cells.size(); i++) {
        msd = cells[i].calc_msd();
        if (cells[i].is_type(H))
        {
            sum_msd_H += msd;
            count_H++;
        }
        else {
            sum_msd_C += msd;
            count_C++;
        }
    }
    fprintf(file, "%.6f\t%.6f\t%.6f\n", time, sum_msd_H/count_H, sum_msd_C/count_C);
}

/// Updates the results file with the current state of the system.
///
/// The output quantities (time, dfc_H, dfc_C) with respect to the scaling
/// quantities ulength and utime.
///
/// Keyword Arguments:
///     file    --  file pointer to the output data location.
///     time    --  current time in units of tau
void Simulation::write_cell_dfc(FILE *file, double time) { 
    double count_H = 0.0, count_C = 0.0;
    double sum_dfc_H = 0.0, sum_dfc_C = 0.0;
    double dfc = 0.0;
    for (unsigned int i=0; i<this->cells.size(); i++) {
        dfc = pow(cells[i].get_rho(),2.0); // square the distance from center
        if (cells[i].is_type(H))
        {
            sum_dfc_H += dfc;
            count_H++;
        }
        else {
            sum_dfc_C += dfc;
            count_C++;
        }
    }
    fprintf(file, "%.6f\t%.6f\t%.6f\n", time, sum_dfc_H/count_H, sum_dfc_C/count_C);
}

/// Finds all cells in the current simulation state that are colliding.
///
/// Current implementation in a brute force O(N^2) direct comparison. Collided 
/// cells are added to the cells' adjacency list for future use in the force
/// calculations. Since this uses the hastable, the method order is really
/// O(N^(3/2))
void Simulation::find_collisions(Table *table) {
    for (int b=0; b<table->num_boxes; b++) { // For each box
        if (table->boxes[b].size() > 1) { // If there are at least two cells
            for (unsigned int i=0; i<table->boxes[b].size(); i++) { // Collide A
                cells[i].adjlst.clear();
                for (unsigned int j=0; j<this->cells.size(); j++) { // with B
                    if (i != j && cells[i].hits(cells[j])) { // If they hit
                            cells[i].adjlst.push_back(cells[j].id);
                    }
                }
            }
        }
    }
}

/// Calculates the adhesive and repulsive forces for each cell.
///
/// Currently, these interactions are set to zero. Future versions will include
/// the JKR model for adhesion and repulsion.
void Simulation::calc_forces() {
    for (unsigned int i=0; i<this->cells.size(); i++) {
        Cell *c1 = &cells[i];
        // Reset Forces
        c1->Fx = 0;
        c1->Fy = 0;
        c1->Fz = 0;
        // Initialize values for JKR force
        double h, E_star, R_star, F, sigma, mag;
        for (int j : c1->adjlst) {
            Cell c2 = cells[j];
            if (c1->is_type(H) && c2.is_type(H)) {
                sigma = sq.surf_E_HH;
            }
            else if (c1->is_type(C) && c2.is_type(C)) {
                sigma = sq.surf_E_CC;
            }
            else {
                sigma = sq.surf_E_HC;
            }
            mag = c1->distance_to(c2);
            h = sq.u_length*(c1->R + c2.R - mag); // put in units of length
            R_star = sq.u_length*(c1->R*c2.R)/(c1->R + c2.R); // put in units of length
            E_star = (4.0/3.0)*(c1->E*c2.E)/
                ((1-pow(c1->nu,2))*c2.E+(1-pow(c2.nu,2))*c1->E);
            // Adh - Rep
            F = sqrt(6*pi*sigma*E_star*pow(R_star,1.5)*pow(h,1.5))
                - E_star*sqrt(R_star)*pow(h,1.5); 
            c1->Fx += F*(c2.x-c1->x)/mag;
            c1->Fy += F*(c2.y-c1->y)/mag;
            c1->Fz += F*(c2.z-c1->z)/mag;
        }
    }
}

/// Update cell positions by integrating the 3D Langevin equations.
///
/// Currently, forces contributing to the particles velocity include Brownian 
/// motion, self-propulsion, cell-cell adhession, and cell-cell repulsion. The
/// numerical integration is done with Forward Euler O(dt).
void Simulation::update_locs() {
    double xnew, ynew, znew;
    for (unsigned int i=0; i<this->cells.size(); i++) {
        #ifdef DEBUG
        if (DEBUG == "brownian") {
            cells[i].dXdt = sqrt(2.0*sq.D/(sq.u_time*dq.dt))*norm_dist(rho_gen);
            cells[i].dYdt = sqrt(2.0*sq.D/(sq.u_time*dq.dt))*norm_dist(rho_gen);
            cells[i].dZdt = sqrt(2.0*sq.D/(sq.u_time*dq.dt))*norm_dist(rho_gen);
        }
        else if (DEBUG == "spp") {
            cells[i].dXdt = (cells[i].self_prop*cos(cells[i].theta)*sin(cells[i].phi))/sqrt(3);
            cells[i].dYdt = (cells[i].self_prop*sin(cells[i].theta)*sin(cells[i].phi))/sqrt(3);
            cells[i].dZdt = (cells[i].self_prop*cos(cells[i].phi))/sqrt(3);
        }
        else { // Debugging just for prints.
            // Apply JKR force 
            cells[i].dXdt = (sq.D/sq.u_energy)*cells[i].Fx;
            cells[i].dYdt = (sq.D/sq.u_energy)*cells[i].Fy;
            cells[i].dZdt = (sq.D/sq.u_energy)*cells[i].Fz;
        
            // Apply self propulsion
            cells[i].dXdt += (cells[i].self_prop*cos(cells[i].theta)*sin(cells[i].phi))/sqrt(3);
            cells[i].dYdt += (cells[i].self_prop*sin(cells[i].theta)*sin(cells[i].phi))/sqrt(3);
            cells[i].dZdt += (cells[i].self_prop*cos(cells[i].phi))/sqrt(3);
        
            // Apply diffusion
            cells[i].dXdt += sqrt(2.0*sq.D/(sq.u_time*dq.dt))*norm_dist(rho_gen);
            cells[i].dYdt += sqrt(2.0*sq.D/(sq.u_time*dq.dt))*norm_dist(rho_gen);
            cells[i].dZdt += sqrt(2.0*sq.D/(sq.u_time*dq.dt))*norm_dist(rho_gen);
        }
        #else // No debug mode specified
        // Apply JKR force 
        cells[i].dXdt = (sq.D/sq.u_energy)*cells[i].Fx;
        cells[i].dYdt = (sq.D/sq.u_energy)*cells[i].Fy;
        cells[i].dZdt = (sq.D/sq.u_energy)*cells[i].Fz;
        
        // Apply self propulsion
        cells[i].dXdt += (cells[i].self_prop*cos(cells[i].theta)*sin(cells[i].phi))/sqrt(3);
        cells[i].dYdt += (cells[i].self_prop*sin(cells[i].theta)*sin(cells[i].phi))/sqrt(3);
        cells[i].dZdt += (cells[i].self_prop*cos(cells[i].phi))/sqrt(3);
        
        // Apply diffusion
        cells[i].dXdt += sqrt(2.0*sq.D/(sq.u_time*dq.dt))*norm_dist(rho_gen);
        cells[i].dYdt += sqrt(2.0*sq.D/(sq.u_time*dq.dt))*norm_dist(rho_gen);
        cells[i].dZdt += sqrt(2.0*sq.D/(sq.u_time*dq.dt))*norm_dist(rho_gen);
        #endif
        
        // Update positions and orientations (forward euler)
        // This is where we remove the dimensions from the force calculations.
        xnew = cells[i].x + dq.dt*(sq.u_time/sq.u_length)*cells[i].dXdt;
        ynew = cells[i].y + dq.dt*(sq.u_time/sq.u_length)*cells[i].dYdt;
        znew = cells[i].z + dq.dt*(sq.u_time/sq.u_length)*cells[i].dZdt;

        // Validity check: stay in sphere
        if ( pow(xnew,2)+pow(ynew,2)+pow(znew,2) < pow(dq.Rb,2))
        {
            cells[i].x = xnew;
            cells[i].y = ynew;
            cells[i].z = znew;
        }

        // Update theta and phi
        cells[i].theta += sqrt(6*sq.D*sq.u_time*dq.dt)*norm_dist(theta_gen)/sq.u_length;
        cells[i].phi += sqrt(6*sq.D*sq.u_time*dq.dt)*norm_dist(theta_gen)/sq.u_length;
    }
}

/// Logic for dividing cells.
/// 
/// Each cell in the simulation is checked for the specified type. If the cell
/// is of the correct type, it is split into two daughter cells of the same 
/// size, each with a velocity half the initial velocity. The two cells are
/// possitioned a distance sigma/2 appart, where sigma is the unit length.
void Simulation::divide_cells(CellType T){
    Cell c, cnew;
    double off = 1.73205; // sqrt(3)*sigma, roughly
    double mag;
    double xn, yn, zn;
    unsigned int N = cells.size();
    for (unsigned int i=0; i<N; i++) {
        if (cells[i].is_type(T)) { // Divide
            c = this->cells[i];
            mag = c.get_rho();
            xn = (1-off/mag)*c.x;
            yn = (1-off/mag)*c.y;
            zn = (1-off/mag)*c.z;
            
            cnew = Cell(pop_size, T, xn, yn, zn, c.theta, c.phi, c.R, c.E, c.nu, 
                    c.self_prop);
            this->cells.push_back(cnew);
            this->pop_size += 1;
        }
    }
}


/// Logic for running the simulation.
/// 
/// Creates the filename used to write out particle position. This includes 
/// opening and closing the file pointer. Order of operations: (1) write out
/// current cell locations, (2) find cells that are colliding, (3) calculate
/// the adhesive and repulsive forces acting on each cell, and (4) update
/// the particle locations.
void Simulation::run() {
    // Get current time
    time_t rawtime;
    struct tm * timeinfo;
    char buffer [80];

    time (&rawtime);
    timeinfo = localtime (&rawtime);
    strftime (buffer,80,"%c",timeinfo);
    
    FILE *loc_file = fopen(this->loc_name, "w");
    FILE *msd_file = fopen(this->msd_name, "w");
    FILE *dfc_file = fopen(this->dfc_name, "w");
    FILE *files[3] = {loc_file, msd_file, dfc_file};

    for (FILE *file : files) {
        fprintf(file, "# Time started: %s\n", buffer); 
        fprintf(file, "# Simulation using %d cells\n", 2*sq.N);
        fprintf(file, "# Bounding Radius: %f\t Time Scale: %f\n", dq.Rb, sq.u_time);
    }
    fprintf(loc_file, "# Format: time(tau)\t cellID\t cellType\t x\t y\t z\n");
    fprintf(msd_file, "# Format: time(tau)\t MSD_H\t MSD_C\n");
    fprintf(dfc_file, "# Format: time(tau)\t DFC_H\t DFC_C\n");

    int total = (int)(dq.tf/dq.dt);
    int data_freq = (int)(total/sq.num_points);
    #ifdef DEBUG
    printf("\n### Starting simulation.\n");
    printf("Taking data every %d steps for %d steps.\n", data_freq, total);
    #endif
    
    // Set up loop counters
    int count = data_freq;      // Used to write out data
    double t = 0;               // Keeps track of time
    double h_div = 0;           // Keeps track of healthy division
    double c_div = 0;           // Keeps track of cancer division
    
    Table table = Table(this);
    while (t < dq.tf) {
        if (count == data_freq) { 
            write_cell_loc(loc_file, t); 
            write_cell_msd(msd_file, t); 
            write_cell_dfc(dfc_file, t); 
            count = 0;
        }
        table.add_cells();
        find_collisions(&table);
        calc_forces();
        update_locs();
        table.clear_table();

        // Check for division
        if (h_div >= dq.hdiv) { 
            printf("HEALTHY DIVISION @ t = %.f\n", t);
            divide_cells(H);
            printf("Current Population size = %d\n", this->pop_size);
            h_div = 0; // reset
        }
        if (c_div >= dq.cdiv) { 
            printf("CANCER  DIVISION @ t = %.f\n", t);
            divide_cells(C);
            printf("Current Population size = %d\n", this->pop_size);
            c_div = 0; // reset
        }

        // Update counters
        t += dq.dt;
        count++;
        h_div += dq.dt;
        c_div += dq.dt;
    }
    write_cell_loc(loc_file, t); 
    write_cell_msd(msd_file, t); 
    write_cell_dfc(dfc_file, t); 
    #ifdef DEBUG
    printf("\n### Simulation Completed.\n");
    #endif

    fclose(loc_file);
    fclose(msd_file);
    fclose(dfc_file);
}
