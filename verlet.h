#ifndef __VERLET_H
#define __VERLET_H

#include "types.h"

void verlet_step1(Positions_t &positions, Velocities_t &velocities,
                  const Forces_t &forces, double timestep);
void verlet_step2(Velocities_t &velocities, const Forces_t &forces, double timestep);

#endif  // __VERLET_H