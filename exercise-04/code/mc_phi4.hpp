#include "4d.hpp"
#include <random>

class phi4 {
protected:
	double kappa_, lambda_;
	std::ranlux48 rng;
	std::uniform_int_distribution<int> randdist_int_0_V;
	std::uniform_real_distribution<double> randdist_double_0_1;

public:
	lattice<double> grid;
	// constructor for L^4 lattice
	phi4 (int L, double kappa, double lambda, unsigned int seed);

	// returns change in energy due to changing phi -> new_phi at site i
	double dE (int i, double new_phi);

	// returns total energy
	double E ();

	// returns potential energy for given phi^2
	double V (double phi2);

	// returns total magnetisation / V
	double M ();

	// does mc update of site, return 1 if update accepted, 0 otherwie
	int mc_update_site (int i);

	// does n_sweeps sweeps, i.e. n_sweeps * VOLUME mc updates,
	// returns acceptance as a double between 0 and 1
	double mc_sweep (int n_sweeps);

	void init_random_field();

	void init_constant_field();

};
