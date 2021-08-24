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


#ifndef YAMD_ATOMS_H
#define YAMD_ATOMS_H

#include "types.h"

class Atoms {
public:
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;
    Mass_t masses;
    Names_t names;
    // number of atoms as input
    Atoms(const int &nb_atoms) :
            positions{3,nb_atoms}, velocities{3,nb_atoms}, forces{3, nb_atoms}, masses{nb_atoms}, names(nb_atoms) {
        positions.setZero();
        velocities.setZero();
        forces.setZero();
        masses.setOnes();
    }
    // positions as input
    Atoms(const Positions_t &p) :
            positions{p}, velocities{3, p.cols()}, forces{3, p.cols()}, masses{p.cols()}, names(p.cols())  {
        velocities.setZero();
        forces.setZero();
        masses.setOnes();
        for(int i=0; i<p.cols(); i++){names[i] = "H";}
    }
    // positions and velocities as input
    Atoms(const Positions_t &p, const Velocities_t &v) :
            positions{p}, velocities{v}, forces{3, p.cols()}, masses{p.cols()}, names(p.cols()) {
        assert(p.cols() == v.cols());
        forces.setZero();
        masses.setOnes();
    }
    // outputs number of atoms
    size_t nb_atoms() const {
        return positions.cols();
    }
};

#endif //YAMD_ATOMS_H
