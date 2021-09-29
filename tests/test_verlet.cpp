/*
* Copyright 2021 Robert Schütze
*
* ### MIT license
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/

#include <gtest/gtest.h>
#include "HeaderFiles/verlet.h"
#include "HeaderFiles/atoms.h"

// test atom movement for simple model problem: constant force and random initial conditions
TEST(VerletTest, BasicAssertions) {
    int n_atoms = 100; // 100 atoms
    Positions_t p = Positions_t::Random(3,n_atoms); // atom positions
    Positions_t p0 = p;
    Velocities_t v = Velocities_t::Random(3,n_atoms); // atom velocities
    Velocities_t v0 = v;
    Forces_t f = Forces_t::Random(3,n_atoms); // atom forces
    double timestep = 0.1;
    int nb_steps = 10000; // 10k time steps
    for (int i = 0; i < nb_steps; ++i) {
        verlet_step1_simple(p, v, f, timestep); // no force update needed
        verlet_step2_simple(v, f, timestep);
    }
    double t = timestep*nb_steps; // total time
    for(int dim = 0; dim<3; dim++){ // check if parabola in all three dimensions
        for(int atom = 0; atom<n_atoms; atom++){
            EXPECT_NEAR(p(dim,atom), p0(dim,atom)+v0(dim,atom)*t+0.5*f(dim,atom)*t*t, 1e-6);
        }
    }
}

// ATOM CLASS edition: test atom movement for simple model problem: constant force and random initial conditions
TEST(VerletTest, BasicAssertionsAtomClass) {
    int n_atoms = 100; // 100 atoms
    Positions_t p = Positions_t::Random(3,n_atoms); // atom positions
    Positions_t p0 = p;
    Velocities_t v = Velocities_t::Random(3,n_atoms); // atom velocities
    Velocities_t v0 = v;
    Forces_t f = Forces_t::Random(3,n_atoms); // atom forces
    Atoms atoms(p0,v0); // use atom class to test corresponding verlet integrator
    atoms.forces = f;
    double timestep = 0.1;
    int nb_steps = 10000; // 10k time steps
    for (int i = 0; i < nb_steps; ++i) {
        verlet_step1(atoms, timestep); // no force update needed
        verlet_step2(atoms, timestep);
    }
    double t = timestep*nb_steps; // total time
    for(int dim = 0; dim<3; dim++){ // check if parabola in all three dimensions
        for(int atom = 0; atom<n_atoms; atom++){
            EXPECT_NEAR(atoms.positions(dim,atom), p0(dim,atom)+v0(dim,atom)*t+0.5*f(dim,atom)*t*t, 1e-6);
        }
    }
}