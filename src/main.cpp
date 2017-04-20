/*
 * File:    main.cpp
 * Author:  Allen Sanford (ras9841@rit.edu)
 * Description:
 *      Deals with setting up, cleaning up, and running the simulation.
 */

// Includes
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <cstring>
#include <sstream>
#include <libconfig.h++>
#include "params.hpp"
#include "simulation.hpp"

/// Prints a message defining proper commandline usage. 
void print_usage() {
    printf("Error with usage: ./run input_file output_file\n"); 
    exit(EXIT_FAILURE);
}

/// Main function for the co-culture simulation
///
/// Deals with reading in the input configuration file, creating the simulation,
/// displaying results to standard out/error, and cleanup. Memory for arrays 
/// containing the input and output filenames is allocated and freed internally.
///
/// Keyword Arguments:
///     argc    --  number of commandline arguments (should be three)
///     argv    --  commandline entries (2nd is input file, 3rd is output file)
int main(int argc, char *argv[]){
    if (argc != 3) { print_usage(); }
	
    std::stringstream in_input;
    std::stringstream in_dfc_out;
    std::stringstream in_msd_out;
    std::stringstream in_loc_out;

    in_input << "../inputs/" << argv[1];
    in_dfc_out << "../outputs/" << argv[2] << "_dfc.csv";
    in_msd_out << "../outputs/" << argv[2] << "_msd.csv";
    in_loc_out << "../outputs/" << argv[2] << "_loc.csv";

    std::string input_str = in_input.str();
    std::string msd_str = in_msd_out.str();
    std::string dfc_str = in_dfc_out.str();
    std::string loc_str = in_loc_out.str();
    
    char *input = new char[input_str.length() + 1];
    char *out_msd = new char[msd_str.length() + 1];
    char *out_dfc = new char[dfc_str.length() + 1];
    char *out_loc = new char[loc_str.length() + 1];

    std::strcpy(input, input_str.c_str());
    std::strcpy(out_msd, msd_str.c_str());
    std::strcpy(out_dfc, dfc_str.c_str());
    std::strcpy(out_loc, loc_str.c_str());

    printf("### Using %s to determine simulation parameters\n", input);
    printf("### Writing cell locations to %s \n", out_loc);
    printf("### Writing cell MSD values to %s \n", out_msd);
    printf("### Writing cell DFC values to %s \n", out_dfc);
    #ifdef DEBUG
    printf("### DEBUG mode: %s\n", DEBUG);
    #endif 

    // (Try to) Generate configuration object from input file.
    libconfig::Config config;
    try {
        config.readFile(input);
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
    Params::SysQuants sq;
    Params::DimensionlessQuants dq;
    double Kb, Temp, pac_frac, r_healthy, r_cancer, surf_E, prop, blk_mod;
    double hdiv, cdiv;
    try {
        // Process scaling parameters
        sq.num_points = config.lookup("Number_Of_Data_Points");
        sq.N = config.lookup("Cells_Per_Population");
        sq.D = config.lookup("Diffusion_Const");
        Kb = config.lookup("Boltzmann_Constant");
        Temp = config.lookup("Temperature");
        sq.u_energy = Kb*Temp; 
        r_healthy = config.lookup("Radius_H");
        r_cancer = config.lookup("Radius_C");
        sq.u_length = r_healthy+r_cancer; // Average diameter
        prop = config.lookup("Propulsion");
        //sq.u_time = sq.u_length*sq.u_length/sq.D;
        sq.u_time = sq.u_length/prop;
        
        // These parameters should have units since the force calculation
        // addresses the scaling.
        sq.radius_H = r_healthy;
        sq.radius_C = r_cancer;
        sq.vp = prop;
        sq.prop_H = prop*(double)config.lookup("Propulsion_H");
        sq.prop_C = prop*(double)config.lookup("Propulsion_C");
        sq.e_mod_H = config.lookup("Elastic_Mod_H");
        sq.e_mod_C = config.lookup("Elastic_Mod_C");
        surf_E = config.lookup("Surface_E");
        sq.surf_E_HH = surf_E*(double)config.lookup("Surface_E_H-H");
        sq.surf_E_HC = surf_E*(double)config.lookup("Surface_E_H-C");
        sq.surf_E_CC = surf_E*(double)config.lookup("Surface_E_C-C");
        blk_mod = config.lookup("Bulk_Mod_H");
        sq.poisson_H = .5-sq.e_mod_H/(6.0*blk_mod);
        blk_mod = config.lookup("Bulk_Mod_C");
        sq.poisson_C = .5-sq.e_mod_C/(6.0*blk_mod);

        // Process remaining settings and create a dimensionless structure
        pac_frac = config.lookup("Packing_Fraction");
        dq.Rb = cbrt(sq.N/(4*pac_frac)); // make unitless by dividing by diam
        dq.dt = config.lookup("Time_Step"); // already unitless
        dq.tf = config.lookup("Simulation_Duration"); // already unitless
        hdiv = (double)config.lookup("Division_Rate_H")*3600;   // now in seconds
        dq.hdiv = hdiv/sq.u_time;                       // now unitless
        cdiv = (double)config.lookup("Division_Rate_C")*3600;   // now in seconds
        dq.cdiv = cdiv/sq.u_time;                       // now unitless
    }
    catch (const libconfig::SettingNotFoundException &snf)
    {
        fprintf(stderr, "Error with configuration file: setting missing.\n");
        exit(EXIT_FAILURE);
    }

    #ifdef DEBUG
    // Print out system params
    printf("### System Parameters\n");
    printf("Cells/Population:\t\t %d\n", sq.N);
    printf("Number of Data Points:\t\t %d\n", sq.num_points);
    printf("Diffusion Const:\t\t %g\n", sq.D);
    printf("Unit Length:\t\t\t %g\n", sq.u_length);
    printf("Unit Energy:\t\t\t %g\n", sq.u_energy);
    printf("Unit Time:\t\t\t %g\n", sq.u_time);
    printf("Healthy Cell Speed:\t\t %g\n", sq.prop_H);
    printf("Cancer Cell Speed:\t\t %g\n", sq.prop_C);
    printf("Healthy Cell Radius:\t\t %g\n", sq.radius_H); 
    printf("Cancer Cell Radius:\t\t %g\n", sq.radius_C); 
    printf("Healthy Cell Elastic Mod:\t %g\n", sq.e_mod_H); 
    printf("Cancer Cell Elastic Mod:\t %g\n", sq.e_mod_C); 
    printf("Healthy Cell Poisson Num:\t %g\n", sq.poisson_H); 
    printf("Cancer Cell Poisson Num:\t %g\n", sq.poisson_C);
    printf("Surface Energy (HH):\t\t %g\n", sq.surf_E_HH); 
    printf("Surface Energy (HC):\t\t %g\n", sq.surf_E_HC); 
    printf("Surface Energy (CC):\t\t %g\n", sq.surf_E_CC); 

    // Print out dimensionless params
    printf("### Dimensionless Parameters\n");
    printf("Packing Fraction:\t\t %g\n", pac_frac);
    printf("Bounding Radius:\t\t %g\n", dq.Rb);
    printf("Time Step:\t\t\t %g\n", dq.dt);    
    printf("Simulation Duration:\t\t %g\n", dq.tf);    
    #endif 

    Simulation sim (out_loc, out_msd, out_dfc, sq, dq);
    sim.run();

    delete [] input;
    delete [] out_msd;
    delete [] out_dfc;
    delete [] out_loc;
    return EXIT_SUCCESS; 
}
  
