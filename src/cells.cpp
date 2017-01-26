#include "cells.hpp"
#include <math.h>

Cell::Cell(int id, CellType type, double x, double y, 
        double z, double thta, double ph){
    this->id = id;
    this->type = type;
    this->x = x;
    this->y = y;
    this->z = z;
    this->theta = thta;
    this->phi = ph;
}

const char *Cell::get_type() {
    return this->type == H ? "H" : "C";
}

int Cell::is_type(CellType tp) {
    return this->type == tp ? 1 : 0;
}

int Cell::hits(Cell c, double diam) {
    if (sqrt(pow(c.x-this->x,2)+pow(c.y-this->y,2)+pow(c.z-this->z,2)) > diam){
        return 0;
    } 
    else {
        return 1;
    }
}
