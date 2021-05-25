//
// Created by schue on 25.05.2021.
//
#include <gtest/gtest.h>
#include <math.h>
#include "verlet.h"
#include "verlet.cpp"

// Demonstrate some basic assertions.
TEST(VerletTest, BasicAssertions) {
    double x = 0.;
    double y = 0.;
    double z = 0.;
    double vx = 0.;
    double vy = 0.;
    double vz = 0.;
    double timestep = 0.1;
    int nb_steps = 100;
    double f = 1.;
    for (int i = 0; i < nb_steps; ++i) {
        verlet_step1(x, y, z, vx, vy, vz, 0., 0., f, timestep);
        verlet_step2(vx, vy, vz, 0., 0., f, 0.1);
    }
    EXPECT_NEAR(z, 0.5*f*pow(timestep*nb_steps,2.), 1e-6);
}
