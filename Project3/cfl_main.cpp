/*	***********************************************************************************************************************************
Computational Optimization in Logistics: Tutorial 1
	1) Reading a Capacitated Facility Location Problem (CFL) instance from file
	2) Storing the data inside an object of a class
	3) Finding a solution by a greedy heuristic
*********************************************************************************************************************************** */


/*	**********************************************************************************************************************************
	Include the file cfldata.h which provides the class for reading and storing the instance data. It also includes other files we need
*********************************************************************************************************************************** */

	#include "cfl_data.h"



/*	***********************************************************************************************************************************
	The main function is the entry point of our program
*********************************************************************************************************************************** */

	int main() {
		cfldata cfl;
		bool quit = !(cfl.read("cfl.dat"));
		if (quit) return 0;
		cflsol sol(cfl);


	}
