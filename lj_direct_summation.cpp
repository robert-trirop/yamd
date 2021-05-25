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

#include "lj_direct_summation.h"
#include "atoms.h"

double lj_direct_summation(Atoms &atoms, double epsilon, double sigma){
    double E = 0.;

    for(int k = 0; k<atoms.nb_atoms(); k++) {
        for (int i = 0; i < atoms.nb_atoms(); i++) {
            if(i != k){
                Eigen::VectorXd rij_v = atoms.positions.col(i)-atoms.positions.col(k);
                double rij = rij_v.norm();
                E += 2.0*epsilon*(pow(sigma/rij,12.)-pow(sigma/rij,6.));
            }

        }
    }
    return E;
}
