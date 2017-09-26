#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "2d.hpp"
#include <complex>

TEST_CASE( "square lattice of ints", "[npoints]" ) {
	lattice<int> grid (16);
	grid.at(1, 2) = 4;
	REQUIRE( grid.at(1,2) == 4 );
	grid.at(-1, 16) = 99;
	REQUIRE( grid.at(15, 0) == 99 );
	grid.at(12, -1) = 43;
	REQUIRE( grid.at(12, 15) == 43 );
	grid.at(-1, -1) = 3;
	REQUIRE( grid.at(15, 15) == 3 );
	grid.at(16, 17) = 7;
	REQUIRE( grid.at(0, 1) == 7 );
}

TEST_CASE( "rectangular lattice of doubles", "[npoints]" ) {
	lattice<double> grid (8, 10);
	grid.at(-1, -1) = 3.14159;
	REQUIRE( grid.at(7, 9) == 3.14159 );
	grid.at(-1, 10) = 2.66;
	REQUIRE( grid.at(7, 0) == 2.66 );
}

//TEST_CASE( "complex double lattice", "[npoints]" ) {
//	lattice<std::complex<double>> grid (4, 4);
//}