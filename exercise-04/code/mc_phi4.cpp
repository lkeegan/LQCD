#include "mc_phi4.hpp"
#include "stats.hpp"
#include <iostream>
#include <vector>
#include <algorithm>

phi4::phi4 (int L, double kappa, double lambda, unsigned int seed) :
		grid (L), kappa_(kappa), lambda_(lambda), rng(seed), 
		randdist_int_0_V (0, L*L*L*L-1), 
		randdist_double_0_1 (0, 1) {
}

double phi4::dE (int i, double new_phi) {
	double phi2 = grid[i]*grid[i];
	double new_phi2 = new_phi*new_phi;
	double nn = (grid.up0(i) + grid.up1(i) + grid.up2(i) + grid.up3(i));
	nn += (grid.dn0(i) + grid.dn1(i) + grid.dn2(i) + grid.dn3(i));
	return -2.0 * kappa_ * (new_phi - grid[i]) * nn + V(new_phi2) - V(phi2);
}

double phi4::V (double phi2) {
	return phi2 * (lambda_ * phi2 + (1.0 - 2.0 * lambda_));
}

double phi4::E () {
	double E = 0;
	for(int i=0; i<grid.V; i++) {
		double phi2 = grid[i]*grid[i];
		E -= 2.0*kappa_*grid[i] * (grid.up0(i) + grid.up1(i) + grid.up2(i) + grid.up3(i));
		E += V(phi2);
	}
	return E / static_cast<double>(grid.V);
}

double phi4::M () {
	double M = 0;
	for(int i=0; i<grid.V; i++) {
		M += grid[i];
	}
	return fabs(M) / static_cast<double>(grid.V);
}

// does mc update of site, return 1 if update accepted, 0 otherwie
int phi4::mc_update_site (int i) {
	double r = randdist_double_0_1 (rng);
	double new_phi = grid[i] + 4.5*(r-0.5);
	double deltaE = static_cast<double>(dE(i, new_phi));
	r = randdist_double_0_1 (rng);
	if(r < exp(-deltaE)) {
		grid[i] = new_phi;
		return 1;
	}
	return 0;
}

// does VOLUME mc updates, returns acceptance as a double between 0 and 1
double phi4::mc_sweep (int n_sweeps) {
	int acc = 0;
	int n_updates = grid.V*n_sweeps;
	for(int i=0; i<n_updates; ++i){
		int j = randdist_int_0_V (rng);
		acc += mc_update_site (j);
	}
	return static_cast<double>(acc) / static_cast<double>(n_updates);
}

// fill field with random distribution in [-0.5,0.5]
void phi4::init_random_field() {
	for(int i=0; i<grid.V; i++) {
		grid[i] = randdist_double_0_1 (rng) - 0.5;
	}
}

// fill field with all +1
void phi4::init_constant_field() {
	for(int i=0; i<grid.V; i++) {
		grid[i] = 1.0;
	}
}

int main(int argc, char *argv[]) {

    if (argc != 6) {
        std::cout << "This program requires 5 arguments:" << std::endl;
        std::cout << "L, kappa, lambda, n_thermalisation, n_measurements" << std::endl;
        std::cout << "e.g. ./MC_phi4 8 0.13 0.02 1e3 1e5" << std::endl;
        return 1;
    }

	int L = atoi(argv[1]);
	double kappa = atof(argv[2]);
	double lambda = atof(argv[3]);
	int n_thermalisation = atof(argv[4]);
	int n_measurements = atof(argv[5]);

	std::cout.precision(17);
	std::cout << "# 4d phi^4 Model MC\t" << std::endl;
	std::cout << "# L\t" << L << std::endl;
	std::cout << "# kappa\t" << kappa << std::endl;
	std::cout << "# lambda\t" << lambda << std::endl;
	std::cout << "# n_thermalisation\t" << n_thermalisation << std::endl;
	std::cout << "# n_measurements\t" << n_measurements << std::endl;

	unsigned int seed = 12345;

	phi4 mc_phi4 (L, kappa, lambda, seed);
	mc_phi4.init_constant_field();

	// thermalisation
	double acc = 0.;
	for(int i=0; i<n_thermalisation; ++i){
		acc += mc_phi4.mc_sweep(1);
	}
	std::cout << "# Acceptance during thermalisation: " << 
		100.0*static_cast<double>(acc)/static_cast<double>(n_thermalisation)
		<< "%" << std::endl;

	// measurements: output average values with standard error
	// note no attempt to measure or account for autocorrelations
	std::vector<double> E, EE, M, MM;
	E.reserve(n_measurements);
	EE.reserve(n_measurements);
	M.reserve(n_measurements);
	MM.reserve(n_measurements);
	acc = 0.;
	for(int i=0; i<n_measurements; ++i){
		acc += mc_phi4.mc_sweep(1);
		double tmpE = mc_phi4.E();
		E.push_back(tmpE);
		EE.push_back(tmpE*tmpE);
		double tmpM = mc_phi4.M();
		M.push_back(tmpM);
		MM.push_back(tmpM*tmpM);
		// optionally output each measurement for later analysis:
		// std::cout << "##" << mc_phi4.E() << "\t" << mc_phi4.M() << std::endl;
	}
	std::cout << "# Acceptance during measurements: " << 
		100.0*static_cast<double>(acc)/static_cast<double>(n_measurements)
		<< "%" << std::endl;
	std::cout << "# kappa" << "\t\tE/V\t\terror" << "\t\t(E/V)^2\t\terror"
			  << "\t\tM/V\t\terror" << "\t\t(M/V)^2\t\terror" << std::endl;
	std::cout << kappa << "\t"
			  << av(E) << "\t" << std_err(E) << "\t"
			  << av(EE) << "\t" << std_err(EE) << "\t"
			  << av(M) << "\t" << std_err(M) << "\t"
			  << av(MM) << "\t" << std_err(MM) << "\t" << std::endl;
	return 0;
}