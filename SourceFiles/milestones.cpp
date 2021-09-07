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
    double mass = 196.96657; // gold atom mass in g/mol or atomic units u (about the same) from http://www.periodensystem.info/elemente/gold/
    // All distances are in Angstrom (10e-10 m)
    double d = 2.885; // atomic distance from reference clusters - corresponds to 408 pm lattice constant
    double cutoff = 10.; // cutoff radius
    double tu = 10.18; // time unit in fs
    // all time variables are multiplied with the unit 10.18 fs to get physical units -> dt = 1.0 -> dt = 10.18 fs
    double dt = 5.0/tu; // timestep
    double t_eq = 100000/tu; // equilibration time
    int steps_eq = t_eq / dt; // number of equilibration steps
    double t_relax = 5000/tu; // relaxation time for each measurement
    int steps_relax = t_relax / dt; // number of relaxation steps
    double dE = 0.002; // energy step in eV per atom
    double tau = 200 * dt; // relaxation time of the thermostat
    int steps = 150; // total number of energy & temperature pair measurements

    std::vector<double> E_list(steps);
    std::vector<double> T_list(steps);

    double E_pot_val;
    double E_kin_val;
    double E_current;
    double E_initial;
    double T_initial;
    double E_target;
    double T_col = 0;
    char filename[50];

    Atoms atoms(0); // initialize atoms object
    int clusters[] = {5, 6, 8, 11, 14, 9};
    for(int c = 0; c<5; c++) {
        // generate cluster atoms object
        atoms = generate_cluster(clusters[c], d, mass);
        std::cout << atoms.nb_atoms() << std::endl;
        NeighborList neighbor_list(cutoff); // initialize neighbor list

        // equilibrate cluster to 500 K
        for (int i = 1; i <= steps_eq; i++) {
            verlet_step1(atoms, dt);
            neighbor_list.update(atoms);
            E_pot_val = gupta(atoms, neighbor_list, cutoff);
            verlet_step2(atoms, dt);
            berendsen_thermostat(atoms, 500, dt, tau);
            //if (i > 3 * steps_eq / 4) { T_col += T(atoms); } // collect temperatures after convergence to check average
            std::cout << i << std::endl;
            if(i%100==0){ // uncomment for trajectory output
                sprintf(filename, "../Data/Out_MS7/traj%04d.xyz", i/100);
                write_xyz(filename, atoms);
            }
        }

        double t_add = 100000/tu; // additional relaxation time without thermostat
        int steps_add = t_add / dt; // additional number of relaxation steps
        for (int i = 1; i <= steps_add; i++) {
            verlet_step1(atoms, dt);
            neighbor_list.update(atoms);
            E_pot_val = gupta(atoms, neighbor_list, cutoff);
            verlet_step2(atoms, dt);
            T_col += T(atoms);
            std::cout << T(atoms) << std::endl;
        }

        E_initial = E_pot_val + E_kin(atoms);
        E_list[0] = E_initial / atoms.nb_atoms();
        T_initial = T_col / steps_add;
        T_list[0] = T_initial;

        // Energy sweep
        for (int n = 1; n < steps; n++) {
            // increase energy by dE
            E_kin_val = E_kin(atoms);
            E_current = E_pot_val + E_kin_val;
            E_target = E_initial + n * dE * atoms.nb_atoms();
            atoms.velocities *= sqrt(1. + (E_target - E_current) / E_kin_val);

            T_col = 0;
            for (int i = 1; i <= steps_relax; i++) {
                verlet_step1(atoms, dt);
                neighbor_list.update(atoms);
                E_pot_val = gupta(atoms, neighbor_list, cutoff);
                verlet_step2(atoms, dt);
                if (i > steps_relax / 5) { T_col += T(atoms); }
            }
            E_list[n] = (E_pot_val + E_kin(atoms)) / atoms.nb_atoms();
            T_list[n] = 5. / 4. * T_col / steps_relax;
            std::cout << n << std::endl;

            sprintf(filename, "../Data/Out_MS7/ES_%04d.xyz", n);
            write_xyz(filename, atoms);
        }

        // Save temperatures and energies
        sprintf(filename, "../Data/Out_MS7/energies_%04d.txt", c);
        std::ofstream e_file(filename);
        sprintf(filename, "../Data/Out_MS7/temperatures_%04d.txt", c);
        std::ofstream t_file(filename);
        e_file << atoms.nb_atoms() << std::endl; // save number of atoms to the energy file
        for (int i = 0; i < steps; ++i) {
            e_file << std::fixed << std::setprecision(32) << E_list[i] << std::endl;
            t_file << std::fixed << std::setprecision(32) << T_list[i] << std::endl;
        }
        e_file.close();
        t_file.close();
    }
    return 0;
}