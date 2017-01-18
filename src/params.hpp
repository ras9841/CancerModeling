#ifndef __PARAMS__
#define __PARAMS__
namespace Params {
    struct ScalingQuants {
        double u_length;    // unit length: cell diameter    
        double u_energy;    // unit energy: Kb*T
        double u_time;      // unit time:   u_length^2/D
        double D;           // diffusion constant
        int N;              // number of cells per population
    };
    struct DimensionlessQuants {
        double Rb;          // radius of the bounding sphere
        double dt;          // time step (tau)
        double tf;          // simulation duration (tau)
        double prop_H;
        double prop_C;
        double rep_H;       
        double rep_C;
        double adh_HH;
        double adh_HC;
        double adh_CC;
    };
}

#endif
