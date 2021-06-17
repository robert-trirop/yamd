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



#include "../HeaderFiles/milestones.h"

// MILESTONE 4
int MS4() {
    double mass = 1.0;
    double epsilon = 1.0;
    double sigma = 1.0;

    double factor = sqrt(mass*sigma*sigma/epsilon);
    double t_tot = 100*factor; // total simulation time
    double dt = 0.01*factor; // timestep (0.03*factor for non constant energy)
    double safe_interval = 1*factor;
    int steps = t_tot/dt; // number of simulation steps
    int safe_index = safe_interval/dt; // safe state every safe_index

    // Read given xyz file and create corresponding Atoms object
    auto [names, positions, velocities]{read_xyz_with_velocities("../Data/lj54.xyz")};
    Atoms atoms(positions, velocities);

    write_xyz("../Data/Out_MS4/traj0000.xyz",atoms);
    double E = lj_direct_summation(atoms, epsilon, sigma) + E_kin(atoms); // initial energy
    std::vector<double> E_list(steps); // list of all energies
    for(int i=1; i<=steps; i++){
        E_list[i-1]=E;
        verlet_step1(atoms, dt);
        E = lj_direct_summation(atoms, epsilon, sigma);
        verlet_step2(atoms, dt);
        E += E_kin(atoms);
        if(i%safe_index == 0){
            char filename[50];
            sprintf(filename, "../Data/Out_MS4/traj%04d.xyz",i/safe_index);
            write_xyz(filename,atoms);
        }
    }
    std::ofstream file("../Data/Out_MS4/energies.txt");
    for(int i=0; i<steps; ++i) {
        file << std::fixed << std::setprecision(32) << E_list[i] << std::endl;
    }
    file.close();
    return 0;
}



// MILESTONE 5
int MS5() {
    double mass = 1.0;
    double epsilon = 1.0;
    double sigma = 1.0;

    double factor = sqrt(mass*sigma*sigma/epsilon);
    double t_tot = 100*factor; // total simulation time
    double dt = 0.01*factor; // timestep (0.03*factor for non constant energy)
    double safe_interval = 1*factor;
    int steps = t_tot/dt; // number of simulation steps
    int safe_index = safe_interval/dt; // safe state every safe_index

    // Read given xyz file and create corresponding Atoms object
    auto [names, positions, velocities]{read_xyz_with_velocities("../Data/lj54.xyz")};
    Atoms atoms(positions, velocities);

    write_xyz("../Data/Out_MS4/traj0000.xyz",atoms);
    double E = lj_direct_summation(atoms, epsilon, sigma) + E_kin(atoms); // initial energy
    std::vector<double> E_list(steps); // list of all energies
    for(int i=1; i<=steps; i++){
        E_list[i-1]=E;
        verlet_step1(atoms, dt);
        E = lj_direct_summation(atoms, epsilon, sigma);
        verlet_step2(atoms, dt);
        E += E_kin(atoms);
        if(i%safe_index == 0){
            char filename[50];
            sprintf(filename, "../Data/Out_MS4/traj%04d.xyz",i/safe_index);
            write_xyz(filename,atoms);
        }
    }
    std::ofstream file("../Data/Out_MS4/energies.txt");
    for(int i=0; i<steps; ++i) {
        file << std::fixed << std::setprecision(32) << E_list[i] << std::endl;
    }
    file.close();
    return 0;
}