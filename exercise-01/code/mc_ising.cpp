#include "mc_ising.hpp"
#include <iostream>
#include <vector>
#include <algorithm>

ising::ising (int L0, int L1, double beta, unsigned int seed) :
		grid (L0, L1), beta_(beta), rng(seed), 
		randdist_int_0_1 (0, 1), 
		randdist_int_0_V (0, L0*L1-1), 
		randdist_double_0_1 (0, 1) {
}

int ising::dE (int i) {
	int h = grid.up0(i);
	h += grid.dn0(i);
	h += grid.up1(i);
	h += grid.dn1(i);
	return 2 * grid[i] * h;
}

int ising::E () {
	int E = 0;
	for(int i=0; i<grid.V; i++) {
		E += grid[i] * (grid.up0(i) + grid.up1(i));
	}
	return -E;
}

int ising::M () {
	int M = 0;
	for(int i=0; i<grid.V; i++) {
		M += grid[i];
	}
	return abs(M);
}

// does mc update of site, return 1 if update accepted, 0 otherwie
int ising::mc_update_site (int i) {
	double deltaE = static_cast<double>(dE(i));
	double r = randdist_double_0_1 (rng);
	if(r < exp(-beta_ * deltaE)) {
		grid.index(i) *= -1;
		return 1;
	}
	return 0;
}

// does VOLUME mc updates, returns acceptance as a double between 0 and 1
double ising::mc_sweep (int n_sweeps) {
	int acc = 0;
	int n_updates = grid.V*n_sweeps;
	for(int i=0; i<n_updates; ++i){
		int j = randdist_int_0_V (rng);
		acc += mc_update_site (j);
	}
	return static_cast<double>(acc) / static_cast<double>(n_updates);
}

// fill field with random distribution of +1 or -1
void ising::init_random_field() {
	for(int i=0; i<grid.V; i++) {
		grid.index(i) = 2 * randdist_int_0_1 (rng) - 1;
	}
}

// fill field with all +1
void ising::init_constant_field() {
	for(int i=0; i<grid.V; i++) {
		grid.index(i) = 1;
	}
}

// return average of values in vector
double av(std::vector<int> &vec) {
	double sum = std::accumulate(vec.begin(), vec.end(), 0.0);
	return sum/static_cast<double>(vec.size());
}

// return std error of average
double std_err(std::vector<int> &vec) {
	double mean = av(vec);

	double s = std::accumulate(vec.begin(), vec.end(), 0.0,
		[mean](double partial_result, double value){ 
			return partial_result + (value - mean) * (value - mean);
		} );
	/*
	// above is equivalent to:
	double s = 0.;
	for (unsigned int i=0; i<vec.size(); i++) {
		s += (vec[i] - mean)*(vec[i] - mean);
	}
	*/
	return sqrt(s)/static_cast<double>(vec.size());
}

int main(int argc, char *argv[]) {

    if (argc != 5) {
        std::cout << "This program requires four arguments:" << std::endl;
        std::cout << "L, beta, n_thermalisation, n_measurements" << std::endl;
        std::cout << "e.g. ./MC_Ising 32 0.44 1e3 1e5" << std::endl;
        return 1;
    }

	int L = atoi(argv[1]);
	double beta = atof(argv[2]);
	int n_thermalisation = atof(argv[3]);
	int n_measurements = atof(argv[4]);

	std::cout.precision(17);
	std::cout << "# Ising Model MC\t" << std::endl;
	std::cout << "# L\t" << L << std::endl;
	std::cout << "# beta\t" << beta << std::endl;
	std::cout << "# n_thermalisation\t" << n_thermalisation << std::endl;
	std::cout << "# n_measurements\t" << n_measurements << std::endl;

	unsigned int seed = 123;

	ising mc_ising (L, L, beta, seed);
	mc_ising.init_random_field();

	// thermalisation
	double acc = 0.;
	for(int i=0; i<n_thermalisation; ++i){
		acc += mc_ising.mc_sweep(1);
	}
	std::cout << "# Acceptance during thermalisation: " << 
		100.0*static_cast<double>(acc)/static_cast<double>(n_thermalisation)
		<< "%" << std::endl;

	// measurements: output average values with standard error
	// note no attempt to measure or account for autocorrelations
	std::vector<int> E, EE, M, MM;
	E.reserve(n_measurements);
	EE.reserve(n_measurements);
	M.reserve(n_measurements);
	MM.reserve(n_measurements);
	acc = 0.;
	for(int i=0; i<n_measurements; ++i){
		acc += mc_ising.mc_sweep(1);
		double tmpE = mc_ising.E();
		E.push_back(tmpE);
		EE.push_back(tmpE*tmpE);
		double tmpM = mc_ising.M();
		M.push_back(tmpM);
		MM.push_back(tmpM*tmpM);
		// optionally output each measurement for later analysis:
		// std::cout << "##" << mc_ising.E() << "\t" << mc_ising.M() << std::endl;
	}
	std::cout << "# Acceptance during measurements: " << 
		100.0*static_cast<double>(acc)/static_cast<double>(n_measurements)
		<< "%" << std::endl;

	double V = static_cast<double>(L*L);
	std::cout << "# beta" << "\t\tE\t\terror" << "\t\tE^2\t\terror"
			  << "\t\tM\t\terror" << "\t\tM^2\t\terror" << std::endl;
	std::cout << beta << "\t"
			  << av(E)/V << "\t" << std_err(E)/V << "\t"
			  << av(EE)/V << "\t" << std_err(EE)/V << "\t"
			  << av(M)/V << "\t" << std_err(M)/V << "\t"
			  << av(MM)/V << "\t" << std_err(MM)/V << "\t" << std::endl;

	return 0;
}