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
*
* (Created by Robert Schütze on 25.05.2021.)
*/

#include "../HeaderFiles/verlet.h"

// Simple version of the Verlet steps not using the Atoms class
void verlet_step1_simple(Positions_t &positions, Velocities_t &velocities,
                  const Forces_t &forces, double timestep){
    velocities += 0.5*timestep*forces;
    positions += timestep*velocities;
}
void verlet_step2_simple(Velocities_t &velocities, const Forces_t &forces, double timestep){
    velocities += 0.5*timestep*forces;
}

// Verlet steps using the Atoms class
void verlet_step1(Atoms &atoms, double timestep){
    atoms.velocities += 0.5*timestep*atoms.forces;
    atoms.positions += timestep*atoms.velocities;
}
void verlet_step2(Atoms &atoms, double timestep){
    atoms.velocities += 0.5*timestep*atoms.forces;
}
