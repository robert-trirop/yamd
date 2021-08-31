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
    double save_interval = 1*factor;
    int steps = t_tot/dt; // number of simulation steps
    int save_index = save_interval/dt; // save state every save_index

    // Read given xyz file and create corresponding Atoms object
    auto [names, positions, velocities]{read_xyz_with_velocities("../Data/lj54.xyz")};
    Atoms atoms(positions, velocities);

    write_xyz("../Data/Out_MS4/traj0000.xyz", atoms);
    double E = lj_direct_summation(atoms, epsilon, sigma) + E_kin(atoms); // initial energy
    std::vector<double> E_list(steps); // list of all energies
    for(int i=1; i<=steps; i++){
        E_list[i-1]=E;
        verlet_step1(atoms, dt);
        E = lj_direct_summation(atoms, epsilon, sigma);
        verlet_step2(atoms, dt);
        E += E_kin(atoms);
        if(i%save_index == 0){
            char filename[50];
            sprintf(filename, "../Data/Out_MS4/traj%04d.xyz", i/save_index);
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
int MS5(bool write_output, int sx, int sy, int sz) {
    double mass = 1.0;
    double epsilon = 1.0;
    double sigma = 1.0;

    double factor = sqrt(mass*sigma*sigma/epsilon);
    double t_tot = 50*factor; // total simulation time
    double dt = 0.02*factor; // timestep (0.03*factor for non constant energy)
    double tau = 10.*dt; // relaxation time
    double save_interval = 0.2*factor; // save every 10 time steps
    int steps = t_tot/dt; // number of simulation steps
    int save_index = save_interval/dt; // save state every save_index

    Atoms atoms(lattice(sx, sy, sz, sigma*1.3));
    if(write_output){
        write_xyz("../Data/Out_MS5/traj0000.xyz", atoms);
    }

    double E = lj_direct_summation(atoms, epsilon, sigma) + E_kin(atoms); // initial energy
    std::vector<double> E_list(steps); // list of all energies

    clock_t start = clock(); // start clock
    for(int i=1; i<=steps; i++){
        E_list[i-1]=E;
        verlet_step1(atoms, dt);
        E = lj_direct_summation(atoms, epsilon, sigma);
        verlet_step2(atoms, dt);
        E += E_kin(atoms);
        berendsen_thermostat(atoms, 200, dt, tau);

        if(i%save_index == 0 && write_output){ // save state every save_index if output active
            char filename[50];
            sprintf(filename, "../Data/Out_MS5/traj%04d.xyz",i/save_index);
            write_xyz(filename, atoms);
            std::cout << T(atoms) << std::endl; // print current temperature
        }
    }
    clock_t stop = clock(); // end clock
    double sim_time = (double) (stop-start)/CLOCKS_PER_SEC;
    std::cout << sim_time << std::endl; // print simulation time

    if(write_output) { // Save energies to txt file
        std::ofstream file("../Data/Out_MS5/energies.txt");
        for (int i = 0; i < steps; ++i) {
            file << std::fixed << std::setprecision(32) << E_list[i] << std::endl;
        }
        file.close();
    }
    return 0;
}


// MILESTONE 6
int MS6(bool write_output, int sx, int sy, int sz) {
    double mass = 1.0;
    double epsilon = 1.0;
    double sigma = 1.0;
    double cutoff = 5.;

    double factor = sqrt(mass*sigma*sigma/epsilon);
    double t_tot = 50*factor; // total simulation time
    double dt = 0.02*factor; // timestep (0.03*factor for non constant energy)
    double tau = 10.*dt; // relaxation time
    double save_interval = factor; // save every 10 time steps
    int steps = t_tot/dt; // number of simulation steps
    int save_index = save_interval/dt; // save state every save_index

    Atoms atoms(lattice(sx, sy, sz, sigma*1.1));
    if(write_output){
        write_xyz("../Data/Out_MS6/traj0000.xyz", atoms);
    }

    NeighborList neighbor_list(cutoff);
    double E = lj(atoms, neighbor_list, epsilon, sigma, cutoff) + E_kin(atoms); // initial energy
    std::vector<double> E_list(steps); // list of all energies

    clock_t start = clock(); // start clock
    for(int i=1; i<=steps; i++){
        E_list[i-1]=E;
        verlet_step1(atoms, dt);
        E = lj(atoms, neighbor_list, epsilon, sigma, cutoff);
        verlet_step2(atoms, dt);
        E += E_kin(atoms);
        berendsen_thermostat(atoms, 200, dt, tau);

        if(i%save_index == 0 && write_output){ // save state every save_index if output active
            char filename[50];
            sprintf(filename, "../Data/Out_MS6/traj%04d.xyz", i/save_index);
            write_xyz(filename, atoms);
            std::cout << T(atoms) << std::endl; // print current temperature
            std::cout << i/save_index << std::endl; // print current file number
        }
    }
    clock_t stop = clock(); // end clock
    double sim_time = (double) (stop-start)/CLOCKS_PER_SEC;
    std::cout << sim_time << std::endl; // print simulation time

    if(write_output) { // Save energies to txt file
        std::ofstream file("../Data/Out_MS6/energies.txt");
        for (int i = 0; i < steps; ++i) {
            file << std::fixed << std::setprecision(32) << E_list[i] << std::endl;
        }
        file.close();
    }
    return 0;
}


// MILESTONE 7
int MS7(bool write_output, int sx, int sy, int sz) {
    double mass = 196.96657; // in g/mol or atomic units u (about the same)
    double cutoff = 5.;

    // all time units are multiplied with 10.18 fs -> dt = 2 -> dt = 20.36 fs
    double dt = 2; // timestep
    double t_eq = 1000; // equilibration time
    int steps_eq = t_eq/dt; // number of equilibration steps
    double t_relax = 50; // relaxation time for each measurement
    int steps_relax = t_relax/dt;
    double tau = 10*dt; // relaxation time of the thermostat

    // load cluster and create atoms object
    auto [names, positions]{read_xyz("../Data/cluster_923.xyz")};
    Atoms atoms(positions, mass);
    NeighborList neighbor_list(cutoff);

    double E_pot_val;
    double E_kin_val;
    double T_col;

    // equilibrate cluster
    for (int i=1; i<=steps_eq; i++) {
        verlet_step1(atoms, dt);
        neighbor_list.update(atoms);
        E_pot_val = gupta(atoms, neighbor_list, cutoff);
        verlet_step2(atoms, dt);
        E_kin_val = E_kin(atoms);
        berendsen_thermostat(atoms, 400, dt, tau);
        if(i>3*steps/4){T_col += T(atoms);} // collect temperatures to check average
    }

    std::cout << T(atoms) << std::endl;

    for (int i=1; i<=steps_relax; i++) {
        verlet_step1(atoms, dt);
        neighbor_list.update(atoms);
        E_pot_val = gupta(atoms, neighbor_list, cutoff);
        verlet_step2(atoms, dt);
        E_kin_val = E_kin(atoms);
    }
    return 0;
}