#ifndef LATTICE_4D_H
#define LATTICE_4D_H

#include<vector>

// 4d lattice with pbcs
// one vector to store data
// another vector stores pointers to nearest neighbours
template<typename T> class lattice {
protected:
	struct nearest_neighbours {
		T* up0;
		T* dn0;
		T* up1;
		T* dn1;
		T* up2;
		T* dn2;
		T* up3;
		T* dn3;
	};
	std::vector<T> data_;
	std::vector<nearest_neighbours> neighbours_;
	int pbcs (int x, int L);
	int index (int x0, int x1, int x2, int x3);
	
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

	// return reference to data at point (x0, x1, x2, x3)
	T& at (int x0, int x1, int x2, int x3);

	// return reference to data with index i
	T& index (int i);

	// also allow [] to return ref to data with index i
	T& operator[](int i) { return index(i); }

	// returns const reference to nearest neighbours
	const T& up0 (int i) { return *neighbours_[i].up0; }
	const T& dn0 (int i) { return *neighbours_[i].dn0; }
	const T& up1 (int i) { return *neighbours_[i].up1; }
	const T& dn1 (int i) { return *neighbours_[i].dn1; }
	const T& up2 (int i) { return *neighbours_[i].up2; }
	const T& dn2 (int i) { return *neighbours_[i].dn2; }
	const T& up3 (int i) { return *neighbours_[i].up3; }
	const T& dn3 (int i) { return *neighbours_[i].dn3; }

	// fill all sites with supplied value
	void fill (T rhs);

};

// constructor
template<typename T> lattice<T>::lattice (int L0, int L1, int L2, int L3) :
 L0(L0), L1(L1), L2(L2), L3(L3) {
	V = L0 * L1 * L2 * L3;
	// Initialise data vector with 1's
	for(int i=0; i<V; ++i){
		data_.push_back(1);
	}
	// Initialise nearest neighbours vector with pointers to neighbours
	data_.reserve(V);
	neighbours_.reserve(V);
	for(int i=0; i<L0; ++i){
		for(int j=0; j<L1; ++j){
			for(int k=0; k<L2; ++k){
				for(int l=0; l<L3; ++l){
					nearest_neighbours tmp_nn;
					tmp_nn.up0 = &data_[index(i+1,j,k,l)];
					tmp_nn.dn0 = &data_[index(i-1,j,k,l)];
					tmp_nn.up1 = &data_[index(i,j+1,k,l)];
					tmp_nn.dn1 = &data_[index(i,j-1,k,l)];
					tmp_nn.up2 = &data_[index(i,j,k+1,l)];
					tmp_nn.dn2 = &data_[index(i,j,k-1,l)];
					tmp_nn.up3 = &data_[index(i,j,k,l+1)];
					tmp_nn.dn3 = &data_[index(i,j,k,l-1)];
					neighbours_.push_back(tmp_nn);
				}
			}
		}
	}
}

// return data at (x0, x1)
template<typename T> T& lattice<T>::at (int x0, int x1, int x2, int x3) {
	return data_[index(x0, x1, x2, x3)];
}

// return data at (x0, x1)
template<typename T> T& lattice<T>::index (int i) {
	return data_[i];
}

// fill all sites with supplied value
template<typename T> void lattice<T>::fill (T rhs) {
	for(int i=0; i<V; i++) {
		data_[i] = rhs;
	}
}

// return global index of point (x0,x1,x2,x3)
template<typename T> int lattice<T>::index(int x0, int x1, int x2, int x3) {
	return pbcs(x0, L0) + L0 * pbcs(x1, L1) + L0 * L1 * pbcs(x2, L2) + L0 * L1 * L2 * pbcs(x3, L3);
}

// apply periodic bcs to point x in [0,L)
template<typename T> int lattice<T>::pbcs(int x, int L) {
	return (x + 2*L) % L;
}
#endif