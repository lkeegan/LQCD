#ifndef LATTICE_1D_H
#define LATTICE_1D_H

#include<vector>

// 1d lattice with pbcs
// one vector to store data
// another vector stores pointers to nearest neighbours
template<typename T> class lattice {
protected:
	struct nearest_neighbours {
		T* up;
		T* dn;
	};
	std::vector<T> data_;
	std::vector<nearest_neighbours> neighbours_;
	int pbcs (int x, int L);
	
public:
	int L;

	// constructor for L-point lattice
	explicit lattice (int L);

	// default destructor is fine:
	//~lattice ();

	// return reference to data with index i
	T& index (int i);

	// also allow [] to return ref to data with index i
	T& operator[](int i) { return index(i); }

	// returns const reference to nearest neighbours
	const T& up (int i) { return *neighbours_[i].up; }
	const T& dn (int i) { return *neighbours_[i].dn; }

	// fill all sites with supplied value
	void fill (T rhs);

};

// constructor
template<typename T> lattice<T>::lattice (int L) : L(L) {
	// Initialise data vector with 1's
	data_.reserve(L);
	for(int i=0; i<L; ++i){
		data_.push_back(1);
	}
	// Initialise nearest neighbours vector with pointers to neighbours
	neighbours_.reserve(L);
	for(int i=0; i<L; ++i){
		nearest_neighbours tmp_nn;
		tmp_nn.up = &data_[pbcs(i+1,L)];
		tmp_nn.dn = &data_[pbcs(i-1,L)];
		neighbours_.push_back(tmp_nn);
	}
}

// return data at site i
template<typename T> T& lattice<T>::index (int i) {
	return data_[i];
}

// fill all sites with supplied value
template<typename T> void lattice<T>::fill (T rhs) {
	for(int i=0; i<L; i++) {
		data_[i] = rhs;
	}
}

// apply periodic bcs to point x in [0,L)
template<typename T> int lattice<T>::pbcs(int x, int L) {
	return (x + 2*L) % L;
}
#endif