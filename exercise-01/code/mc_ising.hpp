#include "2d.hpp"
#include <random>

class ising {
protected:
	lattice<int> grid;
	double beta_;
	std::ranlux48 rng;
	std::uniform_int_distribution<int> randdist_int_0_1;
	std::uniform_int_distribution<int> randdist_int_0_V;
	std::uniform_real_distribution<double> randdist_double_0_1;

public:
	// constructor for L0xL1 lattice
	ising (int L0, int L1, double beta, unsigned int seed);

	// assume L0 = L1 = L if only one length specified
	explicit ising (int L, double beta, unsigned int seed) : 
		ising::ising (L, L, beta, seed) {}

	// returns change in energy due to flipping site (x0,x1)
	int dE (int i);

	// returns total energy
	int E ();

	// returns total magnetisation
	int M ();

	// does mc update of site, return 1 if update accepted, 0 otherwie
	int mc_update_site (int i);

	// does n_sweeps sweeps, i.e. n_sweeps * VOLUME mc updates,
	// returns acceptance as a double between 0 and 1
	double mc_sweep (int n_sweeps);

	void init_random_field();

	void init_constant_field();

};

double H (lattice<int>& grid);