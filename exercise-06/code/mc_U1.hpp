#include "4d.hpp"
#include <random>
#include <complex>
 
// class that contains phases of U[mu], i.e. 4 x double numbers,
// where the corresponding U(1) element is given by exp(i U[mu]) 
class U1 {
private:
	double U_[4];
public:
	// [mu] operator returns U[mu]
	double& operator[](int i) { return U_[i]; }
	double operator[](int i) const { return U_[i]; }
};


class mc_u1 {
protected:
	lattice grid;
	field<U1> U;
	double beta_;
	std::ranlux48 rng;
	std::uniform_int_distribution<int> randdist_int_0_3;
	std::uniform_int_distribution<int> randdist_int_0_V;
	std::uniform_real_distribution<double> randdist_double_0_1;

	// returns change in energy due to changing angle U[ix][mu] from theta to theta_new
	double dE (int ix, int mu, double theta, double theta_new);

	// returns staple of links A around the link U_{\mu}(ix)
	std::complex<double> staple (int ix, int mu);

	// does mc update of link U[ix][mu], return 1 if update accepted, 0 otherwise
	// delta is the max size of the proposed change: delta = [0, 1]
	int mc_update_site (int ix, int mu, double delta);

	// does over-relaxation update of link U[ix][mu]
	void or_update_site (int ix, int mu);

	// returns dirac flux m at ix through mu-nu plane
	// plaquette angle = physical flux in [-pi,pi] + m * 2pi
	int m_mu_nu (int ix, int mu, int nu);

	// returns net dirac flux through mu-nu faces of cube at ix and ix+rho
	int m_flux (int ix, int mu, int nu, int rho);

public:
	// constructor for L^4 lattice
	mc_u1 (int L, double beta, unsigned int seed);

	// returns average plaquette value normalised to [0,1]
	double plaq ();

	// returns average monopole density normalised to [0,1?]
	double monopole ();

	// does n_sweeps MC sweeps, i.e. n_sweeps * VOLUME * 4 mc link updates,
	// returns acceptance as a double between 0 and 1
	double mc_sweep (int n_sweeps, double delta);

	// does n_sweeps OR sweeps, i.e. n_sweeps * VOLUME * 4 over-relaxation link updates,
	void or_sweep (int n_sweeps);

	void init_random_field();

	void init_constant_field();

};