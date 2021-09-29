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

#include "../HeaderFiles/lj_direct_summation.h"

double lj_direct_summation(Atoms &atoms, double epsilon, double sigma){
    double E = 0.;
    atoms.forces.setZero();
    double ssq = sigma*sigma; // only sigma^2 is used (pre-computation)
    for(int j = 0; j<atoms.nb_atoms(); j++) {// Looping over the upper triangle of the pair matrix
        for (int i = j+1; i < atoms.nb_atoms(); i++) { // refer to report for details
                Vector_t rij_v = atoms.positions.col(j)-atoms.positions.col(i);
                double rij_sq = rij_v(0)*rij_v(0)+rij_v(1)*rij_v(1)+rij_v(2)*rij_v(2); // no sqrt, because only rij^2 is needed
                E += 4.*epsilon*(pow(ssq/rij_sq,6.)-pow(ssq/rij_sq,3.)); //no 1/2, because the pairs are not counted twice
                Vector_t Fij = 24.*epsilon*(pow(ssq*rij_sq,3.)-2.*pow(ssq,6.))/pow(rij_sq,7.)*rij_v; //derivation: see report
                //Due to the antisymmetry of the pair force matrix Fij,...
                atoms.forces.col(i) += Fij; // the j-th summand of the i-th force...
                atoms.forces.col(j) -= Fij; // negatively contributes to the j-th force.
        }
    }
    return E;
}


