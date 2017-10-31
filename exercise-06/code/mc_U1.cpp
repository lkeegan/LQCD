#include "mc_U1.hpp"
#include <iostream>
#include <vector>
#include <algorithm>

constexpr double PI = 3.14159265358979323846;
constexpr double TWO_PI = 2.0 * PI;

mc_u1::mc_u1 (int L, double beta, unsigned int seed) :
		grid (L),
		U (grid), 
		beta_(beta), rng(seed), 
		randdist_int_0_3 (0, 3), 
		randdist_int_0_V (0, grid.V-1), 
		randdist_double_0_1 (0, 1) {
}

double mc_u1::dE (int ix, int mu, double theta, double theta_new) {
	std::complex<double> A = staple(ix, mu);
	double theta_staple = std::arg(A);
	double alpha_staple = std::abs(A);
	double deltaE = cos(theta_staple + theta_new) - cos(theta_staple + theta);
	return -beta_ * alpha_staple * deltaE;
}

std::complex<double> mc_u1::staple (int ix, int mu) {
	std::complex<double> A = 0;
	for(int nu=0; nu<4; ++nu) {
		if (nu != mu) {
			int ix_minus_nu = U.idn(ix,nu);
			// up staple
			A += std::polar(1.0, U.up(ix,mu)[nu] - U.up(ix,nu)[mu] - U[ix][nu]);
			// down staple
			A += std::polar(1.0, -U.up(ix_minus_nu,mu)[nu] - U[ix_minus_nu][mu] + U[ix_minus_nu][nu]);
		}
	}
	return A;
}

double mc_u1::plaq () {
	double P = 0;
	for(int ix=0; ix<grid.V; ix++) {
		for(int mu=1; mu<4; ++mu) {
			for(int nu=0; nu<mu; ++nu) {
				P += cos(U[ix][mu] + U.up(ix,mu)[nu] - U.up(ix,nu)[mu] - U[ix][nu]);
			}
		}
	}
	return P / static_cast<double>(6*grid.V);
}

int mc_u1::m_mu_nu (int ix, int mu, int nu) {
	double theta_plaq = U[ix][mu] + U.up(ix,mu)[nu] - U.up(ix,nu)[mu] - U[ix][nu];
	int m = 0;
	while (theta_plaq > PI) {
		++m;
		theta_plaq -= TWO_PI;
	}
	while (theta_plaq < PI) {
		--m;
		theta_plaq += TWO_PI;
	}
	return m;
}

int mc_u1::m_flux (int ix, int mu, int nu, int rho) {
	int ix_plus_rho = U.iup(ix, rho);
	return m_mu_nu (ix_plus_rho, mu, nu) - m_mu_nu (ix, mu, nu);
}

double mc_u1::monopole () {
	double M = 0;
	for(int ix=0; ix<grid.V; ix++) {
		// M_0
		M += fabs(+ m_flux (ix, 1, 2, 3) - m_flux (ix, 1, 3, 2) + m_flux (ix, 2, 3, 1));
		// M_1
		M += fabs(+ m_flux (ix, 0, 2, 3) - m_flux (ix, 0, 3, 2) + m_flux (ix, 2, 3, 0));
		// M_2
		M += fabs(- m_flux (ix, 0, 1, 3) + m_flux (ix, 0, 3, 1) - m_flux (ix, 1, 3, 0));
		// M_3
		M += fabs(+ m_flux (ix, 0, 1, 2) - m_flux (ix, 0, 2, 1) + m_flux (ix, 1, 2, 0));
	}
	return M / static_cast<double>(4*grid.V);
}

// does mc update of gauge link, return 1 if update accepted, 0 otherwie
int mc_u1::mc_update_site (int ix, int mu, double delta) {
	double theta = U[ix][mu];
	double theta_new = theta + TWO_PI * delta * (randdist_double_0_1 (rng) - 0.5);
	// pbs on theta: project into [0,2pi)
	theta_new = fmod((theta_new + 4.0*TWO_PI), TWO_PI);
	double deltaE = dE(ix, mu, theta, theta_new);
	double r = randdist_double_0_1 (rng);
	if(r < exp(-deltaE)) {
		U[ix][mu] = theta_new;
		return 1;
	}
	return 0;
}

// does VOLUME mc updates, returns acceptance as a double between 0 and 1
double mc_u1::mc_sweep (int n_sweeps, double delta) {
	int acc = 0;
	int n_updates = 4*grid.V*n_sweeps;
	for(int i=0; i<n_updates; ++i){
		int ix = randdist_int_0_V (rng);
		int mu = randdist_int_0_3 (rng);
		acc += mc_update_site (ix, mu, delta);
	}
	return static_cast<double>(acc) / static_cast<double>(n_updates);
}

// does or update of gauge link
void mc_u1::or_update_site (int ix, int mu) {
	std::complex<double> A = staple(ix, mu);
	double theta_staple = std::arg(A);
	double theta_new = -U[ix][mu] - 2.0*theta_staple;
	// pbs on theta: project into [0,2pi)
	theta_new = fmod((theta_new + 4.0*TWO_PI), TWO_PI);
	U[ix][mu] = theta_new;
}

// does VOLUME over-relaxation updates
void mc_u1::or_sweep (int n_sweeps) {
	int n_updates = 4*grid.V*n_sweeps;
	for(int i=0; i<n_updates; ++i){
		int ix = randdist_int_0_V (rng);
		int mu = randdist_int_0_3 (rng);
		or_update_site (ix, mu);
	}
}

// fill field with random angles in [,2*pi)
void mc_u1::init_random_field() {
	for(int ix=0; ix<U.V; ix++) {
		for(int mu=0; mu<4; mu++) {
			U[ix][mu] = TWO_PI * randdist_double_0_1 (rng);
		}
	}
}

// fill field with all 0, i.e. unit gauge links
void mc_u1::init_constant_field() {
	for(int ix=0; ix<U.V; ix++) {
		for(int mu=0; mu<4; mu++) {
			U[ix][mu] = 0.0;
		}
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
        std::cout << "L, beta, n_thermalisation, n_measurements, metropolis sweeps / MC update, over-relaxation updates / MC update" << std::endl;
        std::cout << "e.g. ./MC_U1 4 0.1 1e2 1e3 1 5" << std::endl;
        return 1;
    }

	int L = atoi(argv[1]);
	double beta = atof(argv[2]);
	int n_thermalisation = atof(argv[3]);
	int n_measurements = atof(argv[4]);
	int n_metropolis = atof(argv[5]);
	int n_overrelaxation = atof(argv[6]);
	double delta = 0.5; //size of metropolis step: choose beween 0 and 1 to get ~50% acceptance

	std::cout.precision(17);
	std::cout << "# 4d U(1) Monte Carlo\t" << std::endl;
	std::cout << "# L\t" << L << std::endl;
	std::cout << "# beta\t" << beta << std::endl;
	std::cout << "# n_thermalisation\t" << n_thermalisation << std::endl;
	std::cout << "# n_measurements\t" << n_measurements << std::endl;
	std::cout << "# Metropolis sweeps per MC update\t" << n_metropolis << std::endl;
	std::cout << "# Over-relaxation updates per MC update\t" << n_overrelaxation << std::endl;

	unsigned int seed = 123;

	mc_u1 MC_U1 (L, beta, seed);
	MC_U1.init_random_field();

	// thermalisation
	double acc = 0.;
	for(int i=0; i<n_thermalisation; ++i){
		acc += MC_U1.mc_sweep(n_metropolis, delta);
		MC_U1.or_sweep(n_overrelaxation);
	}
	std::cout << "# Acceptance during thermalisation: " << 
		100.0*static_cast<double>(acc)/static_cast<double>(n_thermalisation)
		<< "%" << std::endl;

	// measurements: output average values with standard error
	// some attempt to account for autocorrelations
	std::vector<double> P, PP, M, MM;
	P.reserve(n_measurements);
	PP.reserve(n_measurements);
	M.reserve(n_measurements);
	MM.reserve(n_measurements);
	acc = 0.;
	for(int i=0; i<n_measurements; ++i){
		acc += MC_U1.mc_sweep(n_metropolis, delta);
		MC_U1.or_sweep(n_overrelaxation);

		double tmpP = MC_U1.plaq();
		P.push_back(tmpP);
		PP.push_back(tmpP*tmpP);
		double tmpM = MC_U1.monopole();
		M.push_back(tmpM);
		MM.push_back(tmpM*tmpM);
	}
	std::cout << "# Acceptance during measurements: " << 
		100.0*static_cast<double>(acc)/static_cast<double>(n_measurements)
		<< "%" << std::endl;

	std::cout << "# beta" << "\t\tP\t\terror" << "\t\tP^2\t\terror"
			  << "\t\tM\t\terror" << "\t\tM^2\t\terror" << std::endl;
	std::cout << beta << "\t"
			  << av(P) << "\t" << std_err(P) << "\t"
			  << av(PP) << "\t" << std_err(PP) << "\t"
			  << av(M) << "\t" << std_err(M) << "\t"
			  << av(MM) << "\t" << std_err(MM) << "\t" << std::endl;

	return 0;
}