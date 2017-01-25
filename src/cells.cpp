#include "cells.hpp"
#include <math.h>

Cell::Cell(int id, CellType type, double x, double y, double z){
    this->id = id;
    this->type = type;
    this->x = x;
    this->y = y;
    this->z = z;
    this->vx = 0;
    this->vy = 0;
    this->vz = 0;
}

const char *Cell::get_type() {
    return this->type == H ? "H" : "C";
}
