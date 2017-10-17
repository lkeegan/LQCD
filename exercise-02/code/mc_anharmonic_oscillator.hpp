#include "1d.hpp"
#include <random>

class oscillator {
protected:
	double a_, m0_, mu2_, lambda_, delta_;
	double histogram_width = 6.0;
	std::ranlux48 rng;
	std::uniform_int_distribution<int> randdist_int_0_L;
	std::uniform_real_distribution<double> randdist_double_0_1;

public:
	lattice<double> grid;
	int histogram_n_bins;
	double histogram_bin_width;
	int histogram_measurements = 0;
	std::vector<int> histogram;
	std::vector<double> bins;
	// constructor for L-point lattice
	oscillator (int L, double a, double m0, double mu2, double lambda, unsigned int seed, int nbins);

	// returns change in energy due to setting site i to value x_prime
	double dE (int i, double x_prime);

	// returns total energy
	double E ();

	// returns potential function
	double V (double x);

	// adds current values of x to histogram
	void add_to_histogram ();

	// returns integral of histogram
	double histogram_integral ();

	// returns groundstate energy from virial theorem:
	// E_0 = mu^2 <x^2> + 3 lambda <x^4>
	double E0 ();

	// does mc update of site i, return 1 if update accepted, 0 otherwie
	int mc_update_site (int i);

	// does n_sweeps sweeps, i.e. n_sweeps * VOLUME mc updates,
	// where each mc update is repeated n_hits times per site
	// returns acceptance as a double between 0 and 1
	double mc_sweep (int n_sweeps, int n_hits);

	void init_random_field();

	void init_constant_field (double x);

};

double H (lattice<int>& grid);