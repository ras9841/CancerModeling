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
    double vx;
    double vy;
    double vz;
    double Fx;
    double Fy;
    double Fz;
    std::vector<Cell> adjlst;
    // methods
    Cell(int id, CellType type, double x, double y, double z);
    const char *get_type();
};
#endif
