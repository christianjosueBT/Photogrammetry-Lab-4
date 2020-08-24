// **********************************************************
// **********************************************************
// **********************************************************
// main function
// **********************************************************
// **********************************************************
// **********************************************************


#include "Functions.h"

int main()
{
	// read in ground control points coordinates
	vector<point> coords;
	read("GCPs.txt", coords);

	// read in Image 3314 coordinates
	vector<point> coords_3314;
	string name = "Coords_3314.txt";
	read(name, coords_3314);
	parameters EOP_3314;
	EOP(coords, coords_3314, 3314, EOP_3314);

	// read in image 3316 coordinates
	vector<point> coords_3316;
	name = "Coords_3316.txt";
	read(name, coords_3316);
	parameters EOP_3316;
	EOP(coords, coords_3316, 3316, EOP_3316);

	return 0;
}