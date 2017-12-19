#ifndef LATTICE_SU2_H
#define LATTICE_SU2_H
#include <complex>
#include <Eigen/Dense>
#include <Eigen/StdVector>

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix2cd)
typedef Eigen::Matrix2cd SU2mat;

// class that contains U[mu], i.e. 4 x SU(2) matrices
class gauge {
private:
	SU2mat U_[4];
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	// [mu] operator returns U[mu]
	SU2mat& operator[](int i) { return U_[i]; }
	const SU2mat& operator[](int i) const { return U_[i]; }
};

#endif //LATTICE_SU3_H
