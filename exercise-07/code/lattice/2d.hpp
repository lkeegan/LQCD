#ifndef LATTICE_2D_H
#define LATTICE_2D_H

#include<vector>

// 2d lattice with pbcs
// pair of integer vectors store indices of nearest up/dn neighbours
class lattice {
protected:
	std::vector<int> neighbours_up_;
	std::vector<int> neighbours_dn_;
	int pbcs (int x, int L) const;
	
public:
	int L0;
	int L1;
	int V;
	const int ndim = 2;

	// constructor for L0xL1xL2xL3 lattice
	lattice (int L0, int L1);

	// assume L^2 if only one length specified
	explicit lattice (int L) : lattice::lattice (L, L) {}

	// default destructor is fine:
	//~lattice ();

	// return index to data from 2-vector coordinates (SLOW: for debugging!)
	int index (int x0, int x1) const;

	// returns index of nearest neighbour in mu direction
	int iup (int i, int mu) const { return neighbours_up_[2*i+mu]; }
	int idn (int i, int mu) const { return neighbours_dn_[2*i+mu]; }
};

lattice::lattice (int L0, int L1) :
 L0(L0), L1(L1) {
	V = L0 * L1;
	// Initialise nearest neighbours vector with indices of neighbours
	neighbours_up_.reserve(2*V);
	neighbours_dn_.reserve(2*V);
	for(int j=0; j<L1; ++j){
		for(int i=0; i<L0; ++i){
			neighbours_up_.push_back(index(i+1,j));
			neighbours_up_.push_back(index(i,j+1));
			neighbours_dn_.push_back(index(i-1,j));
			neighbours_dn_.push_back(index(i,j-1));
		}
	}
}

int lattice::index(int x0, int x1) const {
	return pbcs(x0, L0) + L0 * pbcs(x1, L1);
}

int lattice::pbcs(int x, int L) const {
	return (x + 2*L) % L;
}

#endif