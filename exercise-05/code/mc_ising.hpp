#include "4d.hpp"
#include <random>

class ising {
protected:
	lattice grid;
	std::vector<int> spins, cluster_index, cluster_spins, cluster_size;
	double beta_;
	std::ranlux48 rng;
	std::uniform_int_distribution<int> randdist_int_0_1;
	std::uniform_int_distribution<int> randdist_int_0_V;
	std::uniform_real_distribution<double> randdist_double_0_1;

	// returns change in energy due to flipping site (x0,x1)
	int dE (int i);

	// does mc update of site, return 1 if update accepted, 0 otherwie
	int mc_update_site (int i);

	// returns cluster to which spin at site i belongs
	int find_cluster(int i);

	// connect spins at sites x and y, i.e. combine clusters to which they belong
	void connect_clusters (int x, int y);

	// returns maximum cluster size
	int max_cluster_size ();

public:
	// constructor for L^4 lattice
	ising (int L, double beta, unsigned int seed);

	// returns total energy
	int E ();

	// returns total magnetisation
	int M ();

	// does n_sweeps sweeps, i.e. n_sweeps * VOLUME mc updates,
	// returns acceptance as a double between 0 and 1
	double mc_sweep (int n_sweeps);

	// does n_updates cluster updates - returns average of max cluster size / volume for each update 
	double cluster_update (int n_updates);

	void init_random_field();

	void init_constant_field();

};