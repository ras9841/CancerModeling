#include <iostream>
#include <cstdlib>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h> 
#include <libconfig.h++>

#include "params.hpp"

/*! Usage error message. */
void print_usage() {
    printf("Error with usage: ./run filename\n"); 
    exit(EXIT_FAILURE);
}

/*! Main function for the co-culture simulation, */
int main(int argc, char *argv[]){
    if (argc != 2) { print_usage(); }
	
    char input_file[] = "../inputs/"; ///< Sets input directory to ../inputs/ 
    strcat(input_file, argv[1]);
	
    printf("Using %s to determine simulation parameters.\n", input_file);

    // (Try to) Generate configuration object from input file.
    libconfig::Config config;
    try {
        config.readFile(input_file);
    } 
    catch (const libconfig::FileIOException &ferr) {
        fprintf(stderr, "Error with input file IO: invalid read.\n");
        exit(EXIT_FAILURE);
    }
    catch (const libconfig::ParseException &perr) {
        fprintf(stderr, "Error with parsing input: ");
        fprintf(stderr, "(%s : line %d) ", perr.getFile(), perr.getLine());
        fprintf(stderr, "%s\n", perr.getError());
        exit(EXIT_FAILURE);
    }

    // Create parameter structures
    Params::ScalingQuants sq;
    Params::DimensionlessQuants dq;
    double Kb, Temp, pac_frac;
    try {
        // Process scaling parameters
        sq.N = config.lookup("Cells_Per_Population");
        sq.D = config.lookup("Diffusion_Const");
        sq.u_length = config.lookup("Cell_Diameter");
        Kb = config.lookup("Boltzmann_Constant");
        Temp = config.lookup("Temperature");
        sq.u_energy = Kb*Temp; 
        sq.u_time = pow(sq.u_length, 2)/sq.D;

        // Process remaining settings and create a dimensionless structure
        pac_frac = config.lookup("Packing_Fraction");
        dq.Rb = sq.u_length*cbrt(sq.N/(4*pac_frac));
        dq.dt = config.lookup("Time_Step"); // already unitless
        dq.tf = config.lookup("Simulation_Duration"); // already unitless
        dq.prop_H = config.lookup("Propulsion_H");
        dq.prop_C = config.lookup("Propulsion_C");
        dq.rep_H = config.lookup("Repulsion_H");
        dq.rep_C = config.lookup("Repulsion_C");
        dq.adh_HH = config.lookup("Adhesion_H-H");
        dq.adh_HC = config.lookup("Adhesion_H-C");
        dq.adh_CC = config.lookup("Adhesion_C-C");
    }
    catch (const libconfig::SettingNotFoundException &snf)
    {
        fprintf(stderr, "Error with configuration file: setting missing.\n");
        exit(EXIT_FAILURE);
    }

    
    return EXIT_SUCCESS; 
}
  
