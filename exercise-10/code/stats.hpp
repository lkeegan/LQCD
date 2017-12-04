#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>

namespace stats {
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

// return jacknife error of average, ignoring autocorrelations
double jacknife_err0(std::vector<double> &vec) {
	int n = static_cast<int>(vec.size());
	double sum = 0.;
	for (int i=0; i<n; ++i) {
		sum += vec[i];
	}
	std::vector<double> jacknife(n);
	// each jacknife estimate j of the average is
	// the average excluding the j-th data point: 
	for (int i=0; i<n; ++i) {
		jacknife[i] = (sum - vec[i])/static_cast<double>(n-1);
	}
	// and these estimators have the same average as the true average
	// and in the large n limit variance equal to the standard error on this mean:
	double mean = sum / static_cast<double>(n);
	double var = 0.;
	for (int i=0; i<n; ++i) {
		var += (jacknife[i] - mean)*(jacknife[i] - mean);
	}
	var *= static_cast<double>(n-1)/static_cast<double>(n); 
	return sqrt(var);
}

// returns a vector of block-averaged values
// if the original vector has N elements
// the blocked version has M = N/b elements
// where each element is the average of the orgininal elements
// over blocks of size b
// the final block is discarded if it contains less than b elements
std::vector<double> block_average(std::vector<double>& vec, int b) {
	std::vector<double> blocked;
	int N = static_cast<int>(vec.size());
	int M = N/b;
	blocked.reserve(M);
	for (int i_m=0; i_m<M; ++i_m) {
		double av = 0;
		for (int i_b=0; i_b<b; ++i_b) {
			av += vec[i_m*b + i_b];
		}
		av /= static_cast<double>(b);
		blocked.push_back(av);
	}
	return blocked;
}

// returns auto-correlation function Gamma[t], for t up to tmax
// Gamma[t] = \sum_i^{N-t} [ o_{i} o_{i+j} - <o>^2 ]
std::vector<double> auto_corr(std::vector<double>& vec, int tmax) {
	std::vector<double> gamma (tmax);
	double av2 = av(vec);
	av2 = av2 * av2;
	int n = vec.size();

	for (int t=0; t<tmax; ++t) {
		double sum = 0;
		for (int i=0; i<n-t; ++i) {
			sum += vec[i] * vec[i+t] - av2;
		}
		gamma[t] = sum;
	}
	return gamma;
}

// returns running integrated autocorrelation time
std::vector<double> tau_int(std::vector<double>& vec, int tmax) {
	std::vector<double> gamma = auto_corr(vec, tmax);
	std::vector<double> tau(tmax);

	tau[0] = 0.5;
	for (int t=1; t<tmax; ++t) {
		tau[t] = tau[t-1] + gamma[t]/gamma[0];
	}
	return tau;
}

// standard error on average including effect of
// some estimate of autocorrelations in data
// also assigns tau_int and associated error
double std_err(std::vector<double>& vec, double& tau_int, double& tau_int_err) {
	int tmax = vec.size()/50 + 1;
	std::vector<double> gamma = auto_corr(vec, tmax);

	// define tau_int = tau_run[k]
	// for smallest k satisfying k > 6 tau_run[k]
	double tau = 0.5;
	int k = 0;
	while ( (static_cast<double>(k) < 6.0*tau) && (k<tmax) ) {
		++k;
		tau += gamma[k] / gamma[0];
	}
	tau_int = tau;
	tau_int_err = tau * sqrt( static_cast<double>(2*(2*k + 1)) / static_cast<double>(vec.size()) );
	// if we didn't find a k < tmax to satisfy our stability criterion
	// our data series is too short to give a reliable tau_int
	// so set the error on tau_int to infinity to signal this.
	if (k >= tmax-1) {
		tau_int_err = 1.0/0.0;
	}
	double var = 2.0 * tau * gamma[0];
	return(sqrt(var) / static_cast<double>(vec.size()));
}

// standard error on average including effect of
// some estimate of autocorrelations in data
double std_err(std::vector<double>& vec) {
	double ignore0, ignore1;
	return(std_err(vec, ignore0, ignore1));
}

// write average, standard error including tau_int, tau_int to std::cout
void print_av(std::vector<double>& vec, const std::string& name) {
	double tau, tau_err;
	double mean = av(vec);
	double err = std_err(vec, tau, tau_err);
	std::cout << "# " << std::left << std::setw(22) << name
			  << mean << "\t" << err << "\t"
			  << tau << "\t" << tau_err << std::endl;
}

}

