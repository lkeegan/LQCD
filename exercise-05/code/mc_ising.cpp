#include "mc_ising.hpp"
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
	for(int mu=0; mu<4; ++mu) {
		h += spins[grid.iup(i, mu)] + spins[grid.idn(i, mu)];
	}
	return 2 * spins[i] * h;
}

int ising::E () {
	int E = 0;
	for(int i=0; i<grid.V; i++) {
		for(int mu=0; mu<4; ++mu) {
			E += spins[i] * spins[grid.iup(i, mu)];
		}
	}
	return -E;
}

int ising::M () {
	int M = 0;
	for(int i=0; i<grid.V; i++) {
		M += spins[i];
	}
	return abs(M);
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
			for(int mu=0; mu<4; ++mu) {
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

// return average of values in vector
double av(std::vector<double> &vec) {
	double sum = std::accumulate(vec.begin(), vec.end(), 0.0);
	return sum/static_cast<double>(vec.size());
}

// return std error of average, ignoring autocorrelations
double std_err0(std::vector<double> &vec) {
	double mean = av(vec);
	double s = 0.;
	for (unsigned int i=0; i<vec.size(); i++) {
		s += (vec[i] - mean)*(vec[i] - mean);
	}
	return sqrt(s)/static_cast<double>(vec.size());
}

// returns auto-correlation function gamma[t], for t up to tmax
std::vector<double> auto_corr(std::vector<double>& vec, int tmax) {
   std::vector<double> gamma (tmax);

   double sum_lo = 0.0;
   int n = vec.size();
   for (int i=0; i<n; ++i) {
      sum_lo += vec[i];
   }
   double sum_hi = sum_lo;
      
   for (int t=0; t<tmax; ++t) {
      double norm = 1.0/static_cast<double>(n-t);
      double sum = -norm * sum_lo * sum_hi;

      for (int i=t; i<n; i++) {
         sum += vec[i-t] * vec[i];
      }
      
      gamma[t] = norm * sum;
      
      sum_lo -= vec[n-t-1];
      sum_hi -= vec[t];      
   }

   return gamma;
}

// standard error on arverage including effect of (some simple estimate of) autocorrelations in data 
double std_err(std::vector<double>& vec) {
   int tmax = vec.size()/50 + 1;
   std::vector<double> tau(tmax);
   std::vector<double> g = auto_corr(vec, tmax);
   
   tau[0] = 0.5;
   double taumax = 0.0;
   for (int i=1; i<tmax; i++) {
      tau[i] = tau[i-1] + g[i]/g[0];
      if (tau[i]>taumax)
         taumax = tau[i];
      if (tau[i]<=0.0)
         tmax = i;
   }

   double var = 2.0*taumax*g[0]/static_cast<double>(vec.size());
   return(sqrt(var));   
}

int main(int argc, char *argv[]) {

    if (argc != 7) {
        std::cout << "This program requires 6 arguments:" << std::endl;
        std::cout << "L, beta, n_thermalisation, n_measurements, metropolis sweeps / MC update, cluster updates / MC update" << std::endl;
        std::cout << "e.g. ./MC_Ising 32 0.44 1e3 1e5 5 1" << std::endl;
        return 1;
    }

	int L = atoi(argv[1]);
	double beta = atof(argv[2]);
	int n_thermalisation = atof(argv[3]);
	int n_measurements = atof(argv[4]);
	int n_metropolis = atof(argv[5]);
	int n_cluster = atof(argv[6]);

	std::cout.precision(17);
	std::cout << "# 4d Ising Model MC\t" << std::endl;
	std::cout << "# L\t" << L << std::endl;
	std::cout << "# beta\t" << beta << std::endl;
	std::cout << "# n_thermalisation\t" << n_thermalisation << std::endl;
	std::cout << "# n_measurements\t" << n_measurements << std::endl;
	std::cout << "# Metropolis sweeps per MC update\t" << n_metropolis << std::endl;
	std::cout << "# Cluster updates per MC update\t" << n_cluster << std::endl;

	unsigned int seed = 123;

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

	// measurements: output average values with standard error
	// note no attempt to measure or account for autocorrelations
	std::vector<double> E, EE, M, MM;
	E.reserve(n_measurements);
	EE.reserve(n_measurements);
	M.reserve(n_measurements);
	MM.reserve(n_measurements);
	acc = 0.;
	max_cluster_size = 0.;
	double V = static_cast<double>(L*L*L*L);
	for(int i=0; i<n_measurements; ++i){
		acc += mc_ising.mc_sweep(n_metropolis);
		max_cluster_size += mc_ising.cluster_update(n_cluster);
		double tmpE = mc_ising.E()/V;
		E.push_back(tmpE);
		EE.push_back(tmpE*tmpE);
		double tmpM = mc_ising.M()/V;
		M.push_back(tmpM);
		MM.push_back(tmpM*tmpM);
		// optionally output each measurement for later analysis:
		// std::cout << "##" << mc_ising.E() << "\t" << mc_ising.M() << std::endl;
	}
	std::cout << "# Acceptance during measurements: " << 
		100.0*static_cast<double>(acc)/static_cast<double>(n_measurements)
		<< "%" << std::endl;
	std::cout << "# Average of maximum cluster size during thermalisation: " << 
		100.0*max_cluster_size/static_cast<double>(n_measurements)
		<< "%" << std::endl;

	std::cout << "# beta" << "\t\tE\t\terror" << "\t\tE^2\t\terror"
			  << "\t\tM\t\terror" << "\t\tM^2\t\terror" << std::endl;
	std::cout << beta << "\t"
			  << av(E) << "\t" << std_err(E) << "\t"
			  << av(EE) << "\t" << std_err(EE) << "\t"
			  << av(M) << "\t" << std_err(M) << "\t"
			  << av(MM) << "\t" << std_err(MM) << "\t" << std::endl;

	return 0;
}