#include "table.hpp"
#include <math.h>

#define ROUND(a) (int)(0.5+a)

double Table::atan2(double y, double x) {
    //http://math.stackexchange.com/questions/1098487/atan2-faster-approximation
    if (y == 0.0  && x < 0) {
        return this->pi;
    }
    else if (y == 0.0  && x > 0) {
        return 0.0;
    }
    else if (x == 0.0 && y < 0) {
        return 3*this->pi/2;
    }
    else if (x == 0.0 && y > 0) {
        return this->pi/2;
    }
    else if (x == 0.0 && y == 0.0) {
        return this->pi/4;
    }
    else {
        double min_v = abs(x) < abs(y) ? abs(x) : abs(y);
        double max_v = abs(x) > abs(y) ? abs(x) : abs(y);
        double a = min_v/max_v;
        double s = a*a;
        double r = ((-0.0464964749*s+0.15931422)*s-0.327622764)*s*a+a;
        if (abs(y) > abs(x)) {
            return this->pi/2 - r;
        }
        else if (x < 0) {
            return this->pi - r;
        }
        else {
            return -r;
        }
    }
}

double Table::delta_phi(int i) {
    // dt = 1/[(i+1)*sqrt(i^2+2i+2)]
    // N = round(pi/dt)
e   return this->pi/ROUND(this->pi*(i+1.0)*sqrt(i*i+2*i+2.0));
}

double Table::delta_theta(int i) {
    // dt = 1/(i+1)
    // N = round(2pi/dt)
    return 2.0*this->pi/ROUND(2*this->pi*(i+1.0));
}

int Table::num_theta(int i) {
    return ROUND(2*this->pi/delta_theta(i));
}

int Table::num_phi(int i) {
    return ROUND(this->pi/delta_phi(i));
}

Table::Table(Simulation *sim){
    double Rb = sim->dq.Rb;
    this->sim = sim;
    this->scale = 1.0; // one cel diameter
    this->num_rho = (int)(Rb/scale)+1; // shortcut for ceil
    this->num_boxes = 0;
    for (int i=0; i<this->num_rho; i++) {
        this->num_boxes += num_theta(i)*num_phi(i);
    }
    this->boxes = std::vector<std::vector<Cell*>>(num_boxes);
    //this->boxes.reserve(num_boxes);
    //for (int i=0; i<num_boxes; i++) {
    //    boxes[i] = std::vector<Cell*>();
    //}
}

void Table::add_cell(Cell *c) {
    double rho = sqrt(c->x*c->x+c->y*c->y+c->z*c->z);
    double theta = this->atan2(c->y, c->x) + this->pi;
    double phi = acos(c->z/rho);

    int i = (int)(rho/scale);
    int j = (int)(theta/delta_theta(i));
    int k = (int)(phi/delta_phi(i));

    int prev = 0;
    for (int v=0; v<i; v++) {
        prev += num_theta(v)*num_phi(v);
    }

    int index = prev + j*num_phi(i) + k;
    this->boxes[index].push_back(c);
}

void Table::clear_table() {
    for (int i=0; i<num_boxes; i++){
        boxes[i].clear();
    }
}

void Table::add_cells(){
    for (Cell c : sim->cells) {
        add_cell(&c);
    }
}

Table::~Table(){
    for (int i=0; i<num_boxes; i++){
        for (unsigned int j=0; j< boxes[i].size(); j++) {
            delete boxes[i][j];
        }
        boxes[i].clear();
    }
    boxes.clear();
}

