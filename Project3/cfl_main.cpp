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
	#include<limits>

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

		// Make an arbitrary assignment(for each costermer go through facilities
		// (no ordered), assign customer to the first facility if its capacity 
		//is enough for this customer) of customers to the opened facilities
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
		// call the function that tests the feasibility(if all the customers are assigned)
		//and computes sol. cost
		sol.testfeas(cfl);

		// if the solution is feasible, then print it.
		if (sol.feas) {
			std::cout << "Solution feasible!" << std::endl;
			std::cout << "Solution cost: " << sol.cost << std::endl;
			sol.print(cfl, cfl.probname + ".sol");
		}
		else {
			std::cout << "Solution NOT feasible" << std::endl;
		}
		std::cout << "Press Enter to quit" << std::endl;
		//getchar();
		//Until now we can  know that for the tiny_data, the cheapest_facility_subset() 
		//function find a set, that can allow this simple feasible assignment.


		cflsol sol1(cfl);
		gapheu(cfl, OpenPlant, sol1); //cfl is still used, because it doesn't change at all,
			// OpenPlant also must be used, because we need the set that cheapest_facility_subset()
			// function have found to be the initial subset of facilities of this function;
			// but sol can not be used, because it has already contain one solution(one type of 
			// assignment of our problem, so we need create another sol to contain the new 
			// assignment that this function will find.)
		std::cout << "........" << std::endl;
		// call the function that tests the feasibility(if all the customers are assigned)
		//and computes sol. cost
		sol1.testfeas(cfl);

		// if the solution is feasible, then print it.
		if (sol1.feas) {
			std::cout << "Solution feasible!" << std::endl;
			std::cout << "Solution cost: " << sol1.cost << std::endl;
			sol1.print(cfl, cfl.probname + ".sol");
		}
		else {
			std::cout << "Solution NOT feasible" << std::endl;
		}
		std::cout << "Press Enter to quit" << std::endl;
		getchar();
		//Until now we can  know that for the tiny_data, the cheapest_facility_subset() 
		//function find a set, that can allow this simple feasible assignment.
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
		bool feas = true;
		std::cout << "Check initial value of variables needed:" << std::endl; // z solution cost
		std::cout << sol.cost << std::endl; // z solution cost
		std::vector<double> Qbar;
		Qbar.resize(cfl.nPlants);
		for (int j = 0; j < cfl.nPlants; j++) {
			 Qbar[j] = cfl.PlantCapa[j];  //Initial Qbar the rest of Capacity of facilities
			std::cout << Qbar[j] << std::endl;
		}
		
		std::vector<int> U;
		U.resize(cfl.nCustom);
		for (int i = 0; i < cfl.nCustom; i++) {
			U[i] = i;
		}
		
		while (!U.empty() && feas) {

			double delta_star = - std::numeric_limits<double>::infinity();
			int j_star = -1; 
			int i_star = -1; //we will assign customer i_star to facility j_star at the end of this while_loop 

			for (int i = 0; i < U.size(); i++) {

				std::vector<int> F; //the set of index facilities that can serve customer i, temporly not stored
				double delta; //regret of customer i, temporly not stored
				sorter<double> sort;

				for (int j = 0; j < cfl.nPlants; j++) {
					if (OpenPlant[j] && cfl.CustDeman[U[i]] <= Qbar[j]) {
						F.push_back(j);
					}
				}

				if (F.empty()) {
					feas = false;
					break; //this break should jump out of the out_for loop. we don't need to continue checking the rest customers. There exists no solution.
				}
				else if (F.size() == 1) {
					delta = std::numeric_limits<double>::infinity();
				}
				else {
					std::vector<double> CustCosts_F_i;
					for (int j = 0; j < F.size(); j++) {
						CustCosts_F_i.push_back(cfl.CustCosts[U[i]][F[j]]);
					}
					sort.sortperm(CustCosts_F_i, 'l');
					delta = CustCosts_F_i[sort.perm[1]] - CustCosts_F_i[sort.perm[0]];
				}

				//updata delta_star, j_star, i_star
				if (delta > delta_star) {
					delta_star = delta;
					if (F.size() == 1) {
						j_star = F[0];
						i_star = U[i];
					}
					else {
						j_star = F[sort.perm[0]];
						i_star = U[i];
					}
					
				}
			}

			//assign customer i_star to facility j_star
			sol.CusPlant[i_star] = j_star;
			sol.cost += cfl.CustCosts[i_star][j_star];
			Qbar[j_star] -= cfl.CustDeman[i_star];
			for (int i = 0; i < U.size(); i++) {
				if (U[i] == i_star) {
					U[i] = U[U.size()-1];
					break;
				}
			}
			U.resize(U.size() - 1);
			std::cout << i_star << " is assigned to " << j_star << " U size:" << U.size() << std::endl;
		}

	}