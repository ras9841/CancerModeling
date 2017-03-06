/*
 * File:    cells.cpp
 * Author:  Allen Sanford (ras9841@rit.edu)
 * Description:
 *      Contains the class definition for a Cell object and its associated
 *      methods.
 */

// Imports
#include "cells.hpp"
#include <math.h>

/// Constructor for a Cell to be used in the simulation.
/// 
/// Keyword Arguments
///     id      --  unique identification number for the cell
///     type    --  specifies which population the cell is part of
///     x0      --  initial x location of the cell
///     y0      --  initial y location of the cell
///     z0      --  initial z location of the cell
///     theta0  --` initial theta orientation [0,2*pi]
///     phi0    --  initial phi orientation [0,pi]
///     radius
///     e_mod
///     poisson
///     self_prop
Cell::Cell(int id, CellType type, double x0, double y0, 
        double z0, double theta0, double phi0,
	double radius, double e_mod, double poisson,
        double self_prop){
    this->id = id;
    this->type = type;
    this->x0 = x0;
    this->y0 = y0;
    this->z0 = z0;
    this->x = x0;
    this->y = y0;
    this->z = z0;
    this->theta = theta0;
    this->phi = phi0;
    this->R = radius;
    this->E  = e_mod;
    this->nu = poisson;
    this->self_prop = self_prop;
}

/// Returns a string representation of the cell's type.
///
/// Returns
///     H -> healthy cell, C -> cancer cell
const char *Cell::get_type() {
    return this->type == H ? "H" : "C";
}

/// Returns whether the cell is of the specified type.
///
/// Keyword Arguments
///     tp  --  specified cell type
///
/// Returns
///     1 if of the specifed type, 0 otherwise
int Cell::is_type(CellType tp) {
    return this->type == tp ? 1 : 0;
}

/// Determines whether a given cell collides with this cell.
/// 
/// Two cells are said to have collided if the distance between their centers
/// is less than the length scale of the collision.
///
/// Keyword Arguments
///     c       --  cell potentially colliding with this cell.
int Cell::hits(Cell c) {
    double scale = c.R > this->R ? c.R : this->R;
    return distance_to(c) < scale ? 1 : 0;
}

double Cell::calc_msd() {
    return pow(x-x0,2)+pow(y-y0,2)+pow(z-z0,2);
}


double Cell::distance_to(Cell c) {
    return sqrt(pow(x-c.x,2)+pow(y-c.y,2)+pow(z-c.z,2));
}

double Cell::get_rho() {
    return sqrt(pow(x,2)+pow(y,2)+pow(z,2));
}
