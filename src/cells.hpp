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
    int id;
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
    Cell(int id, CellType type, double x, double y, 
            double z, double theta, double phi);
    const char *get_type();
    int hits(Cell c, double diam);
    int is_type(CellType t); 
};
#endif
