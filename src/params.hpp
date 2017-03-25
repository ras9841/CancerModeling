/*
 * File:    params.hpp
 * Author:  Allen Sanford (ras9841@rit.edu)
 * Description:
 *      Defines the system and dimensionless parameters used in the simulaiton.
 *      Both structs are wrapped in the Params namespace for readability.
 */

#ifndef __PARAMS__
#define __PARAMS__
namespace Params {
    struct SysQuants {
        double u_length;    // unit length: average cell diameter    
        double u_energy;    // unit energy: Kb*T
        double u_time;      // unit time:   u_length^2/D
        double D;           // diffusion constant
        int N;              // number of cells per population
        int num_points;     // number of data points taken
        double prop_H;      // cm/s
        double prop_C;      // cm/s
        double radius_H;    // cm
        double radius_C;    // cm
        double e_mod_H;     // g/s^2
        double e_mod_C;     // g/s^2
        double surf_E_HH;   // g/s^2
        double surf_E_HC;   // g/s^2
        double surf_E_CC;   // g/s^2
        double poisson_H;   // 
        double poisson_C;   //
        double vp;          // cm/s
    };
    struct DimensionlessQuants {
        double Rb;          // radius of the bounding sphere
        double dt;          // time step (tau)
        double tf;          // simulation duration (tau)
    };
}

#endif
