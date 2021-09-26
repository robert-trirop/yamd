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
* (Created by Robert Schütze on 09.06.2021.)
*/
#include "../HeaderFiles/useful_functions.h"

// calculate current kinetic energy of the atoms
double E_kin(Atoms &atoms){
    return (atoms.velocities.square().rowwise()*atoms.masses.transpose()).sum()/2;
}

// calculate current temperature based on the kinetic energy of the atoms
double T(Atoms &atoms){
    return 2./3.* E_kin(atoms)/(kBeV*atoms.nb_atoms());
}

double T_LJ(Atoms &atoms){
    return 2./3.* E_kin(atoms)/(atoms.nb_atoms());
}


// inputs are: number of atoms in three dimensions (ny, ny, nz) and the lattice constant sigma
Positions_t lattice(unsigned int nx, unsigned int ny, unsigned int nz, double sigma){
    Positions_t positions(3, nx*ny*nz); // the next 3 lines include all index permutations in 3 dimensions
    positions.row(0) = List_t::LinSpaced(nx,0,nx-1).replicate(ny*nz,1).reshaped(1,nx*ny*nz);
    positions.row(1) = List_t::LinSpaced(ny,0,ny-1).replicate(nz,nx).transpose().reshaped(1,nx*ny*nz);
    positions.row(2) = List_t::LinSpaced(nz,0,nz-1).replicate(1,nx*ny).transpose().reshaped(1,nx*ny*nz);
    return sigma*positions;
}