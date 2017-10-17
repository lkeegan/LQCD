#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "1d.hpp"

TEST_CASE( "1d ints", "[npoints]" ) {
	lattice<int> grid (16);
	grid[0] = 4;
	grid[1] = 2;
	grid[15] = 1;
	REQUIRE( grid.dn(1) == grid[0] );
	REQUIRE( grid.up(15) == grid[0] );
	REQUIRE( grid.up(0) == grid[1] );
	REQUIRE( grid.dn(0) == grid[15] );
}