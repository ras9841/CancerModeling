#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <cstring>
#include <sstream>
#include <libconfig.h++>
#include "params.hpp"
#include "simulation.hpp"

/*! Usage error message. */
void print_usage() {
    printf("Error with usage: ./run input_file output_file\n"); 
    exit(EXIT_FAILURE);
}

/*! Main function for the co-culture simulation, */
int main(int argc, char *argv[]){
    if (argc != 3) { print_usage(); }
	
    std::stringstream in_input;
    std::stringstream in_output;
    in_input << "../inputs/" << argv[1];
    in_output << "../outputs/" << argv[2] << ".csv";
    std::string input_str = in_input.str();
    std::string output_str = in_output.str();
    char *input = new char[input_str.length() + 1];
    char *output = new char[output_str.length() + 1];
    std::strcpy(input, input_str.c_str());
    std::strcpy(output, output_str.c_str());

    printf("### Using %s to determine simulation parameters\n", input);
    printf("### Writing cell locations to %s \n", output);
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
        dq.Rb = sq.u_length*cbrt((4*sq.N)/pac_frac);
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

    #ifdef DEBUG
    // Print out scaling params
    printf("### Scaling Parameters\n");
    printf("Cells/Population:\t %d\n", sq.N);
    printf("Diffusion Const:\t %g\n", sq.D);
    printf("Cell Diameter:\t\t %g\n", sq.u_length);
    printf("Unit Energy:\t\t %g\n", sq.u_energy);
    printf("Unit Time:\t\t %g\n", sq.u_time);

    // Print out dimensionless params
    printf("### Dimensionless Parameters\n");
    printf("Packing Fraction:\t %g\n", pac_frac);
    printf("Bounding Radius:\t %g\n", dq.Rb);
    printf("Time Step:\t\t %g\n", dq.dt);    
    printf("Simulation Duration:\t %g\n", dq.tf);    
    printf("Propulsions [H,C]:\t [%g,%g]\n", dq.prop_H, dq.prop_C);    
    printf("Repulsions [H,C]:\t [%g,%g]\n", dq.rep_H, dq.rep_C);    
    printf("Adhesions [HH,HC,CC]:\t [%g,%g,%g]\n",dq.adh_HH,dq.adh_HC,dq.adh_CC);    
    #endif 

    Simulation sim (output, sq, dq);
    sim.run(sq, dq);

    delete [] input;
    delete [] output;
    return EXIT_SUCCESS; 
}
  
