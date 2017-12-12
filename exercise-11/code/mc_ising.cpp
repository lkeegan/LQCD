#include "mc_ising.hpp"
#include "stats.hpp"
#include <iostream>
#include <vector>
#include <algorithm>

ising::ising (int L, double beta, unsigned int seed) :
		grid (L),
		spins (grid.V), cluster_index (grid.V), cluster_spins (grid.V), cluster_size(grid.V), 
		beta_(beta), rng(seed), 
		randdist_int_0_1 (0, 1), 
		randdist_int_0_V (0, grid.V-1), 
		randdist_double_0_1 (0, 1) {
}

int ising::dE (int i) {
	int h = 0;
	for(int mu=0; mu<grid.ndim; ++mu) {
		h += spins[grid.iup(i, mu)] + spins[grid.idn(i, mu)];
	}
	return 2 * spins[i] * h;;
}

double ising::E () {
	int E = 0;
	for(int i=0; i<grid.V; i++) {
		for(int mu=0; mu<grid.ndim; ++mu) {
			E += spins[i] * spins[grid.iup(i, mu)];
		}
	}
	return -static_cast<double>(E) / static_cast<double>(grid.V);
}

double ising::M () {
	int M = 0;
	for(int i=0; i<grid.V; i++) {
		M += spins[i];
	}
	return static_cast<double>(abs(M)) / static_cast<double>(grid.V);
}

// does mc update of site, return 1 if update accepted, 0 otherwie
int ising::mc_update_site (int i) {
	double deltaE = static_cast<double>(dE(i));
	double r = randdist_double_0_1 (rng);
	if(r < exp(-beta_ * deltaE)) {
		spins[i] *= -1;
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

// returns cluster to which spin at site i belongs
int ising::find_cluster (int i) {
  while (i != cluster_index[i]) {
    // compress path: set parent of i to its grandparent
    cluster_index[i] = cluster_index[cluster_index[i]];
    // climb tree
    i = cluster_index[i];
  }
  return i;
}

// returns largest cluster
int ising::max_cluster_size () {
	int max = 0;
	for(int i=0; i<grid.V; i++) {
		if (cluster_size[i] > max) {
			max = cluster_size[i];
		}
	}
	return max;
}

// connect spins at sites x and y, i.e. combine clusters to which they belong
void ising::connect_clusters (int x, int y) {
  // find current cluster of each spin
  int x_cluster = find_cluster(x);
  int y_cluster = find_cluster(y);

  // if not already connected, then connect them
  if (x_cluster != y_cluster) {
    // connect smaller tree to root of larger tree
    if (cluster_size[x] < cluster_size[y]) {
      cluster_index[x_cluster] = y_cluster;
      cluster_size[y_cluster] += cluster_size[x_cluster];
    } else {
      cluster_index[y_cluster] = x_cluster;
      cluster_size[x_cluster] += cluster_size[y_cluster];
    }
  }
}

// does n_updates cluster updates - returns average cluster size
double ising::cluster_update (int n_updates) {
	double av_max_cluster_size = 0;
	for(int updates=0; updates<n_updates; ++updates) {
		// initially each spin is in its own cluster
		for(int i=0; i<grid.V; i++) {
			cluster_index[i] = i;
			cluster_size[i] = 1;
		}
		// construct clusters
		for(int i=0; i<grid.V; i++) {
			for(int mu=0; mu<grid.ndim; ++mu) {
				int j = grid.iup(i, mu);
				if(spins[i]==spins[j]) {
					double r = randdist_double_0_1 (rng);
					if(r > exp(-2.0 * beta_)) {
 						connect_clusters (i, j);
					}
				}
			}
		}
		// construct random spin for each cluster
		for(int i=0; i<grid.V; i++) {
			cluster_spins[i] = 2 * randdist_int_0_1 (rng) - 1;
		}
		// assign spins to new cluster spin values
		for(int i=0; i<grid.V; i++) {
			spins[i] = cluster_spins[find_cluster(i)];
		}
		av_max_cluster_size += max_cluster_size() / static_cast<double>(n_updates);
	}
	return av_max_cluster_size / static_cast<double>(grid.V);
}

// fill field with random distribution of +1 or -1
void ising::init_random_field() {
	for(int i=0; i<grid.V; i++) {
		spins[i] = 2 * randdist_int_0_1 (rng) - 1;
	}
}

// fill field with all +1
void ising::init_constant_field() {
	for(int i=0; i<grid.V; i++) {
		spins[i] = 1;
	}
}

int main(int argc, char *argv[]) {

    if (argc != 8) {
        std::cout << "This program requires 7 arguments:" << std::endl;
        std::cout << "L, beta, n_thermalisation, n_measurements, metropolis sweeps / MC update, cluster updates / MC update, seed" << std::endl;
        std::cout << "e.g. ./MC_Ising 32 0.44 1e3 1e5 5 1 12345" << std::endl;
        return 1;
    }

	int L = atoi(argv[1]);
	double beta = atof(argv[2]);
	int n_thermalisation = static_cast<int>(atof(argv[3]));
	int n_measurements = static_cast<int>(atof(argv[4]));
	int n_metropolis = atoi(argv[5]);
	int n_cluster = atoi(argv[6]);
	unsigned int seed = atoi(argv[7]);

	std::cout.precision(11);
	std::cout << "# 2d Ising Model MC\t" << std::endl;
	std::cout << "# L\t" << L << std::endl;
	std::cout << "# beta\t" << beta << std::endl;
	std::cout << "# n_thermalisation\t" << n_thermalisation << std::endl;
	std::cout << "# n_measurements\t" << n_measurements << std::endl;
	std::cout << "# Metropolis sweeps per MC update\t" << n_metropolis << std::endl;
	std::cout << "# Cluster updates per MC update\t" << n_cluster << std::endl;


	ising mc_ising (L, beta, seed);
	mc_ising.init_random_field();

	// thermalisation
	double acc = 0.;
	double max_cluster_size = 0.;
	for(int i=0; i<n_thermalisation; ++i){
		acc += mc_ising.mc_sweep(n_metropolis);
		max_cluster_size += mc_ising.cluster_update(n_cluster);
	}
	std::cout << "# Acceptance during thermalisation: " << 
		100.0*static_cast<double>(acc)/static_cast<double>(n_thermalisation)
		<< "%" << std::endl;
	std::cout << "# Average of maximum cluster size during thermalisation: " << 
		100.0*static_cast<double>(max_cluster_size)/static_cast<double>(n_thermalisation)
		<< "%" << std::endl;

	// measurements: output average values with standard error,
	// including effect of autocorrelations
	std::vector<double> E, EE, M, MM;
	E.reserve(n_measurements);
	EE.reserve(n_measurements);
	M.reserve(n_measurements);
	MM.reserve(n_measurements);
	acc = 0.;
	max_cluster_size = 0.;
	for(int i=0; i<n_measurements; ++i){
		// do updates
		acc += mc_ising.mc_sweep(n_metropolis);
		max_cluster_size += mc_ising.cluster_update(n_cluster);
		// do measurements
		double tmpE = mc_ising.E();
		E.push_back(tmpE);
		EE.push_back(tmpE*tmpE);
		double tmpM = mc_ising.M();
		M.push_back(tmpM);
		MM.push_back(tmpM*tmpM);
		// optionally output each measurement for later analysis:
		std::cout << mc_ising.E() << "\t" << mc_ising.M() << std::endl;
	}
	std::cout << "# Acceptance during measurements: " << 
		100.0*static_cast<double>(acc)/static_cast<double>(n_measurements)
		<< "%" << std::endl;
	std::cout << "# Average of maximum cluster size during thermalisation: " << 
		100.0*max_cluster_size/static_cast<double>(n_measurements)
		<< "%" << std::endl;

	/*
	// output running tau_int of E with error estimate:
	double normN = static_cast<double>(E.size());
	std::vector<double> tau = stats::tau_int(E, E.size()/50 + 1);
	std::cout << "# k \t\t tau_run[k] \t\t error" << std::endl;
	for(int i=0; i<static_cast<int>(tau.size()); ++i) {
		double tau_err = tau[i] * sqrt( static_cast<double>(2*(2*i + 1)) / normN );
		std::cout << i << "\t" << tau[i] << "\t" << tau_err << std::endl;
	}
	exit(0);
	*/

	/*
	// output error as a function of block-averaging bin size
	double tauE_, tauE_err_;
	double errorE0 = stats::std_err0(E);
	std::cout << "# error on E without autocorrelations: " << errorE0 << std::endl;
	double errorE = stats::std_err(E, tauE_, tauE_err_);
	std::cout << "# error on E using autocorrelations: " << errorE << "\ttau_E: " << tauE_ << "err: " << tauE_err_ << 	std::endl;
	std::cout << "#block averaged errors on E:" << std::endl;
	for (int b=1; b<static_cast<int>(E.size()/50 + 1); ++b) {
		std::vector<double> blocked_E = stats::block_average(E, b);
		std::cout << b << "\t" << stats::std_err0(blocked_E) << "\t"  << stats::jacknife_err0(blocked_E) << std::endl;
	}
	exit(0);
	*/

	double tauE, tauE_err;
	std::cout << "# beta" << "\t\tE\t\terror" << "\t\tE^2\t\terror"
			  << "\t\tM\t\terror" << "\t\tM^2\t\terror" << "\t\ttau_int(E)\t\terror" << std::endl;
	std::cout << "# " << beta << "\t\t"
			  << stats::av(E) << "\t" << stats::std_err(E, tauE, tauE_err) << "\t"
			  << stats::av(EE) << "\t" << stats::std_err(EE) << "\t"
			  << stats::av(M) << "\t" << stats::std_err(M) << "\t"
			  << stats::av(MM) << "\t" << stats::std_err(MM) << "\t"
			  << tauE << "\t" << tauE_err << "\t" << std::endl;

	return 0;
}