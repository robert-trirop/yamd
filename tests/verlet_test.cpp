//
// Created by Robert Schütze on 25.05.2021.
//
#include <gtest/gtest.h>
#include "verlet.h"
#include "verlet.cpp"

TEST(VerletTest, BasicAssertions) {
    int n_atoms = 100;
    Positions_t p = Positions_t::Random(3,n_atoms);
    Positions_t p0 = p;
    Velocities_t v = Velocities_t::Random(3,n_atoms);
    Velocities_t v0 = v;
    Forces_t f = Forces_t::Random(3,n_atoms);
    double timestep = 0.01;
    int nb_steps = 10000;
     for (int i = 0; i < nb_steps; ++i) {
        verlet_step1(p, v, f, timestep);
        verlet_step2(v, f, timestep);
    }
     double t = timestep*nb_steps;
     for(int dim = 0; dim<3; dim++){
         for(int atom = 0; atom<n_atoms; atom++){
             EXPECT_NEAR(p(dim,atom), p0(dim,atom)+v0(dim,atom)*t+0.5*f(dim,atom)*t*t, 1e-6);
         }
     }
}
