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

/************************************************************************************************************************************
	Function prototype
*************************************************************************************************************************************/
	// this line makes this function visible to the rest of the file.
	void cheapest_facility_subset(cfldata &cfl, std::vector<bool> &OpenPlant); 
	void gapheu(cfldata &cfl, std::vector<bool> &OpenPlant, cflsol &sol);
	
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

		//std::cout << "Press enter to quit." << std::endl;
		//getchar();


		// Belowing code is to test the new function;
		std::vector<bool> OpenPlant(cfl.nPlants, false);
		cheapest_facility_subset(cfl, OpenPlant);

		for (int i = 0; i < cfl.nPlants; i++) {
			std::cout << "opened facilities:" << i << ":"<< OpenPlant[i] << " ";
		}
		//std::cout << "" << std::endl;
		//std::cout << "Press Enter to quit." << std::endl;
		//getchar();

		// Make an arbitrary assignment of customers to the opened facilities
		// Create vector to store the load of each facilities.
		std::vector<double> sumdemands(cfl.nPlants, 0);
		for (int i = 0; i < cfl.nCustom; i++) {
			sol.CusPlant[i] = -1; //initial the customer is currently unassigned.
			for (int j = 0; j < cfl.nPlants; j++) {
				if (OpenPlant[j] && sumdemands[j] + cfl.CustDeman[i] < cfl.PlantCapa[j]) {
					sol.CusPlant[i] = j;
					sumdemands[j] += cfl.CustDeman[i];
					break;
				}
			}
		}
		std::cout << "........" << std::endl;
		// call the function that tests the feasibility and computes sol. cost
		sol.testfeas(cfl);

		// if the solution is feasible, then print it.
		if (sol.feas) {
			std::cout << "Solution feasible!" << std::endl;
			std::cout << "Solution cost" << sol.cost << std::endl;
			sol.print(cfl, cfl.probname + ".sol");
		}
		else {
			std::cout << "Solution NOT feasible" << std::endl;
		}
		std::cout << "Press Enter to quit" << std::endl;
		getchar();
	}
	

	// A function that finds the cheapest subset of facilities that can supply greater or equal
	void cheapest_facility_subset(cfldata &cfl, std::vector<bool> &OpenPlant)  {
		// Here using '&' before the parameters names can ensure that the modifications we make
		// to the input parameters inside the function will be maintained also outside the function;
		// To find the subset of facilities we just use the sorter class from head file to get a permutation 
		// of the facilities that enables to list them for non-decreasing opening cost. Then we set OpenPlan[i]
		// equal to true for first k facilities i = 0,...k such that Q1+...+Qk > Total customer demand;

		// Create a sorter object
		sorter<double> sort;
		sort.sortperm(cfl.OpenCosts, 'l'); //Sort facilities for non-decreasing opening cost;
		double sumdemands = 0;
		double sumcapacity = 0;

		for(int i = 0; i < cfl.nCustom; i++){
			sumdemands += cfl.CustDeman[i];
		}

		for (int j = 0; j < cfl.nPlants; j++) {
			sumcapacity += cfl.PlantCapa[sort.perm[j]];
			OpenPlant[sort.perm[j]] = true;
			if (sumcapacity > sumdemands) {
				break;
			}
		}
	}


	//A function that executes the GAP heuristic using the subset of open facilities 
	//found by cheapest_facility_subset().

	void gapheu(cfldata &cfl, std::vector<bool> &OpenPlant, cflsol &sol) {

	}