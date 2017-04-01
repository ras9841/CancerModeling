/*
 * File:    cells.hpp
 * Author:  Allen Sanford (ras9841@rit.edu)
 * Description:
 *      Lays out the information about the Cell class and defines the 
 *      CellType enum.
 */

#ifndef __CELLS__
#define __CELLS__

#include <vector>

enum CellType { H, C };

class Cell {
    CellType type;
public:
    // Physcal attributes (these have units)
    double R;                   //  radius
    double E;                   //  elastic modulus
    double nu;                  //  poisson number
    double self_prop;           //  surface energy 
    // Other attributes
    int id;
    double x0;
    double y0;
    double z0;
    double x;
    double y;
    double z;
    double Fx;
    double Fy;
    double Fz;
    double dXdt;
    double dYdt;
    double dZdt;
    double theta;   // [0,2*pi] (vp orientation)
    double phi;     // [0,pi]   (vp orientation)
    std::vector<int> adjlst;    // List of adjacent cells (collided)
    // methods
    Cell() {};
    Cell(int id, CellType type, double x, double y, 
            double z, double theta, double phi,
	    double radius, double e_mod, double poisson,
	    double self_prop);
    const char *get_type();
    int hits(Cell c);
    int is_type(CellType t); 
    double distance_to(Cell c);
    double get_rho();
    double calc_msd();
};
#endif
