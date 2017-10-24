#ifndef LATTICE_4D_H
#define LATTICE_4D_H

#include<vector>

// 4d lattice with pbcs
// pair of integer vectors store indices of nearest up/dn neighbours
class lattice {
protected:
	std::vector<int> neighbours_up_;
	std::vector<int> neighbours_dn_;
	int pbcs (int x, int L) const;
	
public:
	int L0;
	int L1;
	int L2;
	int L3;
	int V;

	// constructor for L0xL1xL2xL3 lattice
	lattice (int L0, int L1, int L2, int L3);

	// assume L^4 if only one length specified
	explicit lattice (int L) : lattice::lattice (L, L, L, L) {}

	// default destructor is fine:
	//~lattice ();

	// return index to data from 4-vector coordinates (SLOW: for debugging!)
	int index (int x0, int x1, int x2, int x3) const;

	// returns index of nearest neighbour in mu direction
	int iup (int i, int mu) const { return neighbours_up_[4*i+mu]; }
	int idn (int i, int mu) const { return neighbours_dn_[4*i+mu]; }
};

lattice::lattice (int L0, int L1, int L2, int L3) :
 L0(L0), L1(L1), L2(L2), L3(L3) {
	V = L0 * L1 * L2 * L3;
	// Initialise nearest neighbours vector with indices of neighbours
	neighbours_up_.reserve(4*V);
	neighbours_dn_.reserve(4*V);
	for(int l=0; l<L3; ++l){
		for(int k=0; k<L2; ++k){
			for(int j=0; j<L1; ++j){
				for(int i=0; i<L0; ++i){
					neighbours_up_.push_back(index(i+1,j,k,l));
					neighbours_up_.push_back(index(i,j+1,k,l));
					neighbours_up_.push_back(index(i,j,k+1,l));
					neighbours_up_.push_back(index(i,j,k,l+1));
					neighbours_dn_.push_back(index(i-1,j,k,l));
					neighbours_dn_.push_back(index(i,j-1,k,l));
					neighbours_dn_.push_back(index(i,j,k-1,l));
					neighbours_dn_.push_back(index(i,j,k,l-1));
				}
			}
		}
	}
}

int lattice::index(int x0, int x1, int x2, int x3) const {
	return pbcs(x0, L0) + L0 * pbcs(x1, L1) + L0 * L1 * pbcs(x2, L2) + L0 * L1 * L2 * pbcs(x3, L3);
}

int lattice::pbcs(int x, int L) const {
	return (x + 2*L) % L;
}

#endif