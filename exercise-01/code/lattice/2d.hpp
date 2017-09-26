#ifndef LATTICE_2D_H
#define LATTICE_2D_H

#include<vector>

// 2d lattice with pbcs
// one vector to store data
// another vector stores pointers to nearest neighbours
template<typename T> class lattice {
protected:
	struct nearest_neighbours {
		T* up0;
		T* dn0;
		T* up1;
		T* dn1;
	};
	std::vector<T> data_;
	std::vector<nearest_neighbours> neighbours_;
	int pbcs (int x, int L);
	int index (int x0, int x1);
	
public:
	int L0;
	int L1;
	int V;

	// constructor for L0xL1 lattice
	lattice (int L0, int L1);

	// assume L0 = L1 = L if only one length specified
	explicit lattice (int L) : lattice::lattice (L, L) {}

	// default destructor is fine:
	//~lattice ();

	// return reference to data at point (x0, x1)
	T& at (int x0, int x1);

	// return reference to data with index i
	T& index (int i);

	// also allow [] to return ref to data with index i
	T& operator[](int i) { return index(i); }

	// returns const reference to nearest neighbours
	const T& up0 (int i) { return *neighbours_[i].up0; }
	const T& dn0 (int i) { return *neighbours_[i].dn0; }
	const T& up1 (int i) { return *neighbours_[i].up1; }
	const T& dn1 (int i) { return *neighbours_[i].dn1; }

	// fill all sites with supplied value
	void fill (T rhs);

};

// constructor
template<typename T> lattice<T>::lattice (int L0, int L1) :
 L0(L0), L1(L1) {
	V = L0 * L1;
	// Initialise data vector with 1's
	for(int i=0; i<V; ++i){
		data_.push_back(1);
	}
	// Initialise nearest neighbours vector with pointers to neighbours
	data_.reserve(V);
	neighbours_.reserve(V);
	for(int i=0; i<L0; ++i){
		for(int j=0; j<L1; ++j){
			nearest_neighbours tmp_nn;
			tmp_nn.up0 = &data_[index(i+1,j)];
			tmp_nn.dn0 = &data_[index(i-1,j)];
			tmp_nn.up1 = &data_[index(i,j+1)];
			tmp_nn.dn1 = &data_[index(i,j-1)];
			neighbours_.push_back(tmp_nn);
		}
	}
}

// return data at (x0, x1)
template<typename T> T& lattice<T>::at (int x0, int x1) {
	return data_[index(x0, x1)];
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

// return global index of point (x0,x1)
template<typename T> int lattice<T>::index(int x0, int x1) {
	return pbcs(x0, L0) + L0 * pbcs(x1, L1);
}

// apply periodic bcs to point x in [0,L)
template<typename T> int lattice<T>::pbcs(int x, int L) {
	return (x + 2*L) % L;
}
#endif