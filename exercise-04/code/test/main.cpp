#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "4d.hpp"
#include <complex>

TEST_CASE( "square lattice of ints", "[npoints]" ) {
	lattice<int> grid (16);
	grid.at(1, 2, 0, 0) = 4;
	REQUIRE( grid.at(1, 2, 0, 0) == 4 );
	grid.at(-1, 16, 5, 7) = 99;
	REQUIRE( grid.at(15, 0, 5, 7) == 99 );
	grid.at(2, 3, 16, -1) = 43;
	REQUIRE( grid.at(2, 3, 0, 15) == 43 );
	grid.at(-1, -1, -1, -1) = 3;
	REQUIRE( grid.at(15, 15, 15, 15) == 3 );
	grid.at(16, 17, 16, 17) = 7;
	REQUIRE( grid.at(0, 1, 0, 1) == 7 );
}

TEST_CASE( "8x10x3x5 lattice of doubles", "[npoints]" ) {
	lattice<double> grid (8, 10, 3, 5);
	grid.at(-1, -1, -1, -1) = 3.14159;
	REQUIRE( grid.at(7, 9, 2, 4) == 3.14159 );
	grid.at(-1, 10, 4, 5) = 2.66;
	REQUIRE( grid.at(7, 0, 1, 0) == 2.66 );
}
