#include "4d.hpp"
#include "su2.hpp"
#include <random>
#include <complex>
 
class mc_su2 {
protected:
	lattice grid;
	field<gauge> U;
	double beta_;
	std::ranlux48 rng;
	std::uniform_int_distribution<int> randdist_int_0_3;
	std::uniform_int_distribution<int> randdist_int_0_V;
	std::uniform_real_distribution<double> randdist_double_0_1;

	// returns staple of links A around the link U_{\mu}(ix)
	SU2mat staple (int ix, int mu);

	// does heat bath update of link U[ix][mu]
	void hb_update_site (int ix, int mu);

	// does over-relaxation update of link U[ix][mu]
	void or_update_site (int ix, int mu);

	// converts real 4-vector of a's into SU(2) matrix
	SU2mat U_from_a (double a0, double a1, double a2, double a3);

public:
	// constructor for T x L^3 lattice
	mc_su2 (int T, int L, double beta, unsigned int seed);

	// returns average plaquette value normalised to 1
	double plaq ();

	// returns average polyakov loop in time direction normalised to 1
	double polyakov_loop ();

	// returns average LxL wilson loop normalised to 1
	double wilson_loop (int L);

	// does n_sweeps HB sweeps, i.e. n_sweeps * VOLUME * 4 heat-bath link updates,
	void hb_sweep (int n_sweeps);

	// does n_sweeps OR sweeps, i.e. n_sweeps * VOLUME * 4 over-relaxation link updates,
	void or_sweep (int n_sweeps);

	void init_random_field ();

	void init_unit_field ();
};