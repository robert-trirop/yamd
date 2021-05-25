#include "verlet.h"

void verlet_step1(double &x, double &y, double &z, double &vx, double &vy, double &vz,
                  double fx, double fy, double fz, double timestep) {
    vx += timestep*0.5*fx;
    vy += timestep*0.5*fy;
    vz += timestep*0.5*fz;
    x  += timestep*vx;
    y  += timestep*vy;
    z  += timestep*vz;
}

void verlet_step2(double &vx, double &vy, double &vz, double fx, double fy, double fz,
                  double timestep) {
    vx += timestep*0.5*fx;
    vy += timestep*0.5*fy;
    vz += timestep*0.5*fz;
}