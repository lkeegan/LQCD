#include "mc_anharmonic_oscillator.hpp"
#include <iostream>
#include <vector>
#include <algorithm>

oscillator::oscillator (int L, double a, double m0, double mu2, double lambda, unsigned int seed, int nbins=50) :
		grid (L), a_(a), m0_(m0), mu2_(mu2), lambda_(lambda), rng(seed), 
		randdist_int_0_L (0, L-1), 
		randdist_double_0_1 (0, 1),
		histogram_n_bins(nbins),
		histogram (nbins, 0), bins(nbins) {
			delta_ = 2.0 * sqrt(a_);
			// bins is the central value of the bin
			// takes values in range [-histogram_width, +histogram_width]
			for(int i=0; i<nbins; ++i) {
				bins[i] = 2.0*histogram_width*((i+0.5)/static_cast<double>(nbins)-0.5);
			}
}

void oscillator::add_to_histogram() {
	for(int i=0; i<grid.L; i++) {
		int bin = static_cast<int>(static_cast<double>(histogram_n_bins)*(grid[i]/(2.0*histogram_width)+0.5));
		//std::cout << grid[i] << "\t" << bin << std::endl;
		++histogram[bin];
	}
	++histogram_measurements;
}

double oscillator::histogram_integral () {
	double bin_width = 2.0*histogram_width/static_cast<double>(histogram_n_bins);
	return bin_width * static_cast<double>(histogram_measurements * grid.L);
}

double oscillator::dE (int i, double x_prime) {
	double x = grid[i];
	double x_updn = grid.up(i) + grid.dn(i); 
	double delE = x_prime*(x_prime - x_updn) - x*(x - x_updn);
	return (m0_ / a_) * delE + a_* (V(x_prime) - V(x));
}

double oscillator::E () {
	double E = 0;
	for(int i=0; i<grid.L; i++) {
		double x = grid[i];
		E += 0.5 * (m0_ / a_) * pow(grid.up(i) - x, 2) + a_ * V(x);
	}
	return E;
}

double oscillator::E0 () {
	double E0 = 0;
	for(int i=0; i<grid.L; i++) {
		double x = grid[i];
		double x2 = x*x;
		E0 += mu2_*x2 + 3.0*lambda_*x2*x2;
	}
	return E0;
}

double oscillator::V (double x) {
	double x2 = x*x;
	return 0.5*mu2_*x2 + lambda_*x2*x2;
}

// does mc update of site, return 1 if update accepted, 0 otherwie
int oscillator::mc_update_site (int i) {
	double r = randdist_double_0_1 (rng);
	double x_prime = grid[i] + 2.0*(r-0.5)*delta_;
	double deltaE = dE(i, x_prime);
	r = randdist_double_0_1 (rng);
	if(r < exp(-deltaE)) {
		grid[i] = x_prime;
		return 1;
	}
	return 0;
}

// does VOLUME mc updates, returns acceptance as a double between 0 and 1
double oscillator::mc_sweep (int n_sweeps, int n_hits) {
	int acc = 0;
	int n_updates = grid.L*n_sweeps;
	for(int i=0; i<n_updates; ++i){
		int site = randdist_int_0_L (rng);
		for(int j=0; j<n_hits; ++j){
			acc += mc_update_site (site);
		}
	}
	return static_cast<double>(acc) / static_cast<double>(n_updates * n_hits);
}

// fill field with flat random values in [0,1)
void oscillator::init_random_field() {
	for(int i=0; i<grid.L; i++) {
		grid[i] = randdist_double_0_1 (rng);
	}
}

// fill field with constant value
void oscillator::init_constant_field (double x) {
	for(int i=0; i<grid.L; i++) {
		grid[i] = x;
	}
}

// return average of values in vector
double av(std::vector<double> &vec) {
	double sum = std::accumulate(vec.begin(), vec.end(), 0.0);
	return sum/static_cast<double>(vec.size());
}

// return std error of average
double std_err(std::vector<double> &vec) {
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

    if (argc != 7) {
        std::cout << "This program requires 6 arguments:" << std::endl;
        std::cout << "L, mu, lambda, n_thermalisation, n_measurements, seed" << std::endl;
        std::cout << "e.g. ./MC_Oscillator 128 0.8 -1.2 1e3 1e5 12345" << std::endl;
        return 1;
    }

	int L = atoi(argv[1]);
	double mu2 = atof(argv[2]);
	double lambda = atof(argv[3]);
	int n_thermalisation = atof(argv[4]);
	int n_measurements = atof(argv[5]);
	double a = 0.1;
	double m0 = 1.0;
	unsigned int seed = atoi(argv[6]);

	std::cout.precision(17);
	std::cout << "# Ising Model MC\t" << std::endl;
	std::cout << "# L\t" << L << std::endl;
	std::cout << "# a\t" << a << std::endl;
	std::cout << "# m0\t" << m0 << std::endl;
	std::cout << "# mu^2\t" << mu2 << std::endl;
	std::cout << "# lambda\t" << lambda << std::endl;
	std::cout << "# n_thermalisation\t" << n_thermalisation << std::endl;
	std::cout << "# n_measurements\t" << n_measurements << std::endl;

	oscillator mc_oscillator (L, a, m0, mu2, lambda, seed);
	mc_oscillator.init_random_field();

	// thermalisation
	double acc = 0;
	for(int i=0; i<n_thermalisation; ++i){
		acc += mc_oscillator.mc_sweep(5, 5);
	}
	std::cout << "# Acceptance during thermalisation: " << 
		100.0*static_cast<double>(acc)/static_cast<double>(n_thermalisation)
		<< "%" << std::endl;

	// measurements: output average values with standard error
	// note no attempt to measure or account for autocorrelations
	std::vector<double> E0;
	std::vector<std::vector<double>> corr(L/2);
	E0.reserve(n_measurements);
	for(int j=0; j<L/2; j++){
		corr[j].reserve(n_measurements);
	}
	acc = 0.;
	for(int i=0; i<n_measurements; ++i){
		acc += mc_oscillator.mc_sweep(5, 5);
		mc_oscillator.add_to_histogram();
		E0.push_back(mc_oscillator.E0());
		if(i%10 == 0) {
			for(int j=0; j<L/2; j++){
				corr[j].push_back(mc_oscillator.grid[j]*mc_oscillator.grid[0]);
			}
		}
	}
	std::cout << "# Acceptance during measurements: " << 
		100.0*static_cast<double>(acc)/static_cast<double>(n_measurements)
		<< "%" << std::endl;

	// output histogram
	for(int i=0; i<mc_oscillator.histogram_n_bins; ++i){
		std::cout << "#hist " << mc_oscillator.bins[i] << "\t" 
			<< static_cast<double>(mc_oscillator.histogram[i])/mc_oscillator.histogram_integral() << std::endl;
	}

	// output E_0
	double V = static_cast<double>(L);
	std::cout << "# mu2" << "\t\tE_0\t\terror" << std::endl;
	std::cout << "#" << mu2 << "\t" << av(E0)/V << "\t" << std_err(E0)/V << std::endl;

	// output correlator
	for(int j=0; j<L/2; j++){
		std::cout << j << "\t" << av(corr[j]) << "\t" << std_err(corr[j]) << std::endl;
	}

	return 0;
}