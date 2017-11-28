#include "mc_SU2.hpp"
#include "stats.hpp"
#include <iostream>
#include <vector>
#include <algorithm>

constexpr double PI = 3.14159265358979323846;
constexpr double TWO_PI = 2.0 * PI;

mc_su2::mc_su2 (int L, double beta, unsigned int seed) :
		grid (L),
		U (grid), 
		beta_(beta), rng(seed), 
		randdist_int_0_3 (0, 3), 
		randdist_int_0_V (0, grid.V-1), 
		randdist_double_0_1 (0, 1) {
}

SU2mat mc_su2::staple (int ix, int mu) {
	SU2mat A = SU2mat::Zero();
	int ix_plus_mu = U.iup(ix,mu);
	for(int nu=0; nu<4; ++nu) {
		if(nu!=mu) {
			A += U[ix_plus_mu][nu] * U.up(ix,nu)[mu].adjoint() * U[ix][nu].adjoint();
			A += U.dn(ix_plus_mu,nu)[nu].adjoint() * U.dn(ix,nu)[mu].adjoint() * U.dn(ix,nu)[nu];
		}
	}
	return A;
}

double mc_su2::plaq () {
	double P = 0;
	for(int ix=0; ix<grid.V; ix++) {
		for(int mu=1; mu<4; ++mu) {
			for(int nu=0; nu<mu; ++nu) {
				P += ((U[ix][mu]*U.up(ix,mu)[nu])*((U[ix][nu]*U.up(ix,nu)[mu]).adjoint())).trace().real();
			}
		}
	}
	return P / static_cast<double>(2*6*grid.V);
}

double mc_su2::polyakov_loop () {
	double p = 0;
	for(int ix3=0; ix3<U.VOL3; ++ix3) {
		int ix = U.it_ix(0, ix3);
		SU2mat P = U[ix][0];
		for(int x0=1; x0<U.L0; x0++) {
			ix = U.iup(ix, 0);
			P *= U[ix][0];
		}
		p += P.trace().real();
	}
	return p / static_cast<double>(2 * U.VOL3);
}


double mc_su2::wilson_loop (int L) {
	double W = 0;
	for(int ix=0; ix<grid.V; ix++) {
		for(int mu=1; mu<4; ++mu) {
			for(int nu=0; nu<mu; ++nu) {
				SU2mat V = SU2mat::Identity();
				int ix_L = ix;
				for(int rho=0; rho<L; ++rho) {
					V *= U[ix_L][mu];
					ix_L = U.iup(ix_L, mu);
				}
				for(int rho=0; rho<L; ++rho) {
					V *= U[ix_L][nu];
					ix_L = U.iup(ix_L, nu);
				}
				for(int rho=0; rho<L; ++rho) {
					ix_L = U.idn(ix_L, mu);
					V *= (U[ix_L][mu]).adjoint();
				}
				for(int rho=0; rho<L; ++rho) {
					ix_L = U.idn(ix_L, nu);
					V *= (U[ix_L][nu]).adjoint();
				}
				W += V.trace().real();
			}
		}
	}
	return W / static_cast<double>(2*6*grid.V);
}

SU2mat mc_su2::U_from_a (double a0, double a1, double a2, double a3) {
	SU2mat A = SU2mat::Zero();
	A(0,0) = std::complex<double>(a0, a3);
	A(0,1) = std::complex<double>(a2, a1);
	A(1,0) = std::complex<double>(-a2, a1);
	A(1,1) = std::complex<double>(a0, -a3);
	return A;
}

void mc_su2::hb_update_site (int ix, int mu) {
	SU2mat A = staple(ix, mu);
	double k = sqrt(A.determinant().real());
	// construct a0 from distribution p(a0) ~ sqrt(1-a0^2) exp(beta k a0)
	std::uniform_real_distribution<double> randdist_z (exp(-beta_ * k), exp(beta_ * k));
	double z;
	double a0;
	double r;
	do {
		// first sample from distribution p(a0) ~ exp(beta k a0)
		z = randdist_z (rng);
		a0 = log(z) / (beta_ * k);
		r = sqrt(1.0 - a0*a0);
	}
	// then accept with probability sqrt(1-a0^2) to give desired distribution
	while (r < randdist_double_0_1 (rng));
	// construct a1, a2, a3 for 3-vector of length r pointing in a random direction
	// random theta in [0, pi]
	double cos_theta = randdist_double_0_1 (rng) * 2.0 - 1.0;
	double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
	// random phi in [0, 2pi]
	double phi = TWO_PI * randdist_double_0_1 (rng);
	double a1 = r * sin_theta * cos(phi);
	double a2 = r * sin_theta * sin(phi);
	double a3 = r * cos_theta;
	// convert a's to SU2 matrix
	SU2mat V = U_from_a (a0, a1, a2, a3);
	U[ix][mu] = (1.0/k) * V * A.adjoint();
}

void mc_su2::hb_sweep (int n_sweeps) {
	int n_updates = 4*grid.V*n_sweeps;
	for(int i=0; i<n_updates; ++i){
		int ix = randdist_int_0_V (rng);
		int mu = randdist_int_0_3 (rng);
		hb_update_site (ix, mu);
	}
}

void mc_su2::or_update_site (int ix, int mu) {
	SU2mat A = staple(ix, mu);
	double k2 = A.determinant().real();
	U[ix][mu] = (1.0/k2)*(A*U[ix][mu]*A).adjoint();
}

void mc_su2::or_sweep (int n_sweeps) {
	int n_updates = 4*grid.V*n_sweeps;
	for(int i=0; i<n_updates; ++i){
		int ix = randdist_int_0_V (rng);
		int mu = randdist_int_0_3 (rng);
		or_update_site (ix, mu);
	}
}

void mc_su2::init_random_field() {
	for(int ix=0; ix<U.V; ix++) {
		for(int mu=0; mu<4; mu++) {
			// make random U(2) matrix from 4 random a's in range [-1,1]
			double a0 = -1.0 + 2.0*randdist_double_0_1 (rng);
			double a1 = -1.0 + 2.0*randdist_double_0_1 (rng);
			double a2 = -1.0 + 2.0*randdist_double_0_1 (rng);
			double a3 = -1.0 + 2.0*randdist_double_0_1 (rng);
			SU2mat A = U_from_a(a0, a1, a2, a3);
			// normalise determinant to 1, i.e. project into SU(2)
			double k = sqrt(A.determinant().real());
			U[ix][mu] = A / k;
		}
	}
}

void mc_su2::init_unit_field() {
	for(int ix=0; ix<U.V; ix++) {
		for(int mu=0; mu<4; mu++) {
			U[ix][mu] = SU2mat::Identity();
		}
	}
}

int main(int argc, char *argv[]) {

    if (argc != 7) {
        std::cout << "This program requires 6 arguments:" << std::endl;
        std::cout << "L, beta, n_thermalisation, n_measurements, heatbath sweeps / MC update, over-relaxation updates / MC update" << std::endl;
        std::cout << "e.g. ./MC_SU2 4 0.1 1e2 1e3 1 5" << std::endl;
        return 1;
    }

	int L = atoi(argv[1]);
	double beta = atof(argv[2]);
	int n_thermalisation = atof(argv[3]);
	int n_measurements = atof(argv[4]);
	int n_heatbath = atof(argv[5]);
	int n_overrelaxation = atof(argv[6]);

	std::cout.precision(10);
	std::cout << "# 4d SU(2) Monte Carlo\t" << std::endl;
	std::cout << "# L\t" << L << std::endl;
	std::cout << "# beta\t" << beta << std::endl;
	std::cout << "# n_thermalisation\t" << n_thermalisation << std::endl;
	std::cout << "# n_measurements\t" << n_measurements << std::endl;
	std::cout << "# Heatbath sweeps per MC update\t" << n_heatbath << std::endl;
	std::cout << "# Over-relaxation updates per MC update\t" << n_overrelaxation << std::endl;

	unsigned int seed = 123;

	mc_su2 MC_SU2 (L, beta, seed);
	MC_SU2.init_random_field();
	//MC_SU2.init_unit_field();

	// thermalisation
	std::cout << "# plaquette \t 4x4 wilson loop" << std::endl;
	for(int i=0; i<n_thermalisation; ++i){
		std::cout << MC_SU2.plaq() << "\t" << MC_SU2.wilson_loop(4) << std::endl;
		MC_SU2.hb_sweep(n_heatbath);
		MC_SU2.or_sweep(n_overrelaxation);
	}

	// measurements: output all measurements for later analysis
	std::cout << "# polyakov_loop, 1x1, 2x2, 3x3, 4x4, 5x5, 6x6" << std::endl;
	for(int i=0; i<n_measurements; ++i){
		MC_SU2.hb_sweep(n_heatbath);
		MC_SU2.or_sweep(n_overrelaxation);
		std::cout << MC_SU2.polyakov_loop() << "\t"
				  << MC_SU2.wilson_loop(1) << "\t"
				  << MC_SU2.wilson_loop(2) << "\t"
				  << MC_SU2.wilson_loop(3) << "\t"
				  << MC_SU2.wilson_loop(4) << "\t"
				  << MC_SU2.wilson_loop(5) << "\t"
				  << MC_SU2.wilson_loop(6) << "\t" << std::endl;
	}

	return 0;
}