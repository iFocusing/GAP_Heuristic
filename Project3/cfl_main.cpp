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

		std::cout << "Capacitated Facility Location Problem Name:" << cfl.probname << std::endl;
		std::cout << "The number of facilities:" << cfl.nPlants << std::endl;
		std::cout << "The number of customers:" << cfl.nCustom << std::endl; 

		double SumOpenCosts = 0;
		double SumPlantCapa = 0;
		double SumCustCosts = 0;

		for (int j = 0; j < cfl.nPlants; j++) {
			SumOpenCosts = SumOpenCosts + cfl.OpenCosts[j];
			SumPlantCapa = SumPlantCapa + cfl.PlantCapa[j];
			for (int i = 0; i < cfl.nCustom; i++){
				SumCustCosts += cfl.CustCosts[i][j];
			}
		}

		std::cout << "Average open cost:" << SumOpenCosts / cfl.nPlants << std::endl;
		std::cout << "Average plant capacity:" << SumPlantCapa / cfl.nPlants << std::endl;
		std::cout << "Average custome cost:" << SumCustCosts / cfl.nCustom << std::endl;

		std::cout << "Press enter to quit." << std::endl;
		getchar();

	}
