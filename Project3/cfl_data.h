
/*	***********************************************************************************************************************************
	Include files: Provide access to functions classes constants and macros defined in include files filename.h :
	use <filename> to include a file located in one of the locations specified by the compiler 
	use "filename.h" to include a file located in the same directory of this file 
	*********************************************************************************************************************************** */

	#include <fstream>			// for file manipulation
	#include <iostream>			// for console input output
	#include <string>			// for strings manipulation
	#include <vector>			// for using vectors (arrays with several built in functionalities)
	#include <algorithm>		// for using various algorithms incl. sort 

	#pragma once

/*	***********************************************************************************************************************************
	Define a class to read and store the instance data. We assume that the input file is formatted as follows:
	-------------------------------------------------
	name_of_instance  #_of_facilities  #_of_customers
	demand_of_customer_i (i=1,...,#_of_customers)
	capacity_of_facility_j (j=1,...,#_of_facilities)
	opening_cost_of_facility_j (j=1,...,#_of_facilities)
	supply_cost_from_facility_1_to_customer_i (i=1,...,#_of_customers)
	supply_cost_from_facility_2_to_customer_i (i=1,...,#_of_customers)
	...
	-------------------------------------------------
	*********************************************************************************************************************************** */

	class cfldata {
	public:
		int nPlants;								// # of facilities (a.k.a plants)
		int nCustom;								// # of customers
		std::vector<double> OpenCosts;				// Facility opening costs 
		std::vector<double> CustDeman;				// Customers demands
		std::vector<double> PlantCapa;				// Facility capacities	
		std::vector<std::vector<double>> CustCosts;	// Customer service costs
		std::string probname;						// A name for the instance

	//  Constructor function: must have same name as the class
	//	Is called automatically to initialize a cfldata object upon creation
		cfldata() {
			nPlants = 0;
			nCustom = 0;
			probname = "";
		}

	//  read: Reads a CFL instance from file. 
	//	In input a string with the location (path) of the instance file
		bool read(std::string path) {	
		
		//  Create an input stream attached to the input file
			std::ifstream input(path);
			if (input.is_open()) {

			//  Read instance name, # of facilities, # of customers
				input >> probname;
				input >> nPlants >> nCustom;

			//  Dimension the vectors 
				OpenCosts.resize(nPlants);
				CustDeman.resize(nCustom);
				PlantCapa.resize(nPlants);
				CustCosts.resize(nCustom, std::vector<double>(nPlants));

			//  Read customer demands
				for (int i = 0; i < nCustom; ++i) {
					input >> CustDeman[i];
					CustDeman[i] = std::floor(CustDeman[i]);
				}
			//  Read facility capacities
				for (int j = 0; j < nPlants; ++j) {
					input >> PlantCapa[j];
					PlantCapa[j] = std::floor(PlantCapa[j]);
				}
			//  Read customer facility open costs
				for (int j = 0; j < nPlants; ++j) {
					input >> OpenCosts[j];
				}
			//  Read customer service costs
				for (int j = 0; j < nPlants; ++j) {
					for (int i = 0; i < nCustom; i++) {
						input >> CustCosts[i][j];
					}
				}
			//  Close the input file
				input.close();
			}
			else { return false; }

			return true;
		}
	};



/*	***********************************************************************************************************************************
	Define a class to store print and check feasibility of a cfl solution
	*********************************************************************************************************************************** */

	class cflsol {
	public:
		double cost;					// Solution cost
		bool feas;						// True if solution is feasible
		std::vector<int> CusPlant;		// CusPlant[i]: The facility customer i is assigned to 

	//  Constructor function: must have same name as the class
	//	Is called automatically to initialize a cflsol object upon creation
		cflsol(cfldata &cfl) {
			cost = 0;
			CusPlant.resize(cfl.nCustom,-1);
		}

	//  testfeas:  checks if the solution stored is feasible
	//  and compute its cost. Gets in input a cfldata object  
	//	returns true if solution feasible or false otherwise
		bool testfeas(cfldata &cfl) {
			std::vector<double> PlantLoad(cfl.nPlants,0);
			std::vector<bool> PlantOpen(cfl.nPlants, false);

			feas = true;
			cost = 0;

    //		Run over each customer, and compute the load of the facilities
			for (int i = 0; i < cfl.nCustom; i++) {
				int j = CusPlant[i];
				if (j < 0 || j >= cfl.nPlants) {
			//		customer i is not assigned
					cost = std::numeric_limits<double>::infinity();
					feas = false;
					break;
				}
				else {
					PlantLoad[j] += cfl.CustDeman[i];
					PlantOpen[j] = true;
					cost += cfl.CustCosts[i][j];
					if (PlantLoad[j] > cfl.PlantCapa[j]) {
				//		Facility j is overloaded
						cost = std::numeric_limits<double>::infinity();
						feas = false;
						break;
					}
				}
			}
			for (int j = 0; j < cfl.nPlants; j++) {
				if (PlantOpen[j]) cost += cfl.OpenCosts[j];
			}
			return feas;
		}

	//  printsol: prints the solution stored to a file named "filename"
		void print(cfldata &cfl, std::string filename) {
			std::ofstream sol(filename);
			if (sol.is_open()) {
				sol << "CFL instance " << cfl.probname << std::endl;
				sol << "Cost " << cost << std::endl << "Feasible " << feas << std::endl;
				sol << "Customers assigned to each facility: " << std::endl << std::endl;
				for (int j = 0; j < cfl.nPlants; j++) {
					bool first_of_facility = true;					
					for (int i = 0; i < cfl.nCustom; i++) {
						if (CusPlant[i] == j) {
							if (first_of_facility) {
								sol << "Facility " << j << std::endl;
								first_of_facility = false;
							}
							sol << i << " ";
						}
					}
					if(!first_of_facility) sol << std::endl << std::endl;
				}
				sol.close();
			}
		}
	};



/*	***********************************************************************************************************************************
	Define a class to sort a vector of values. Values are of any Type
	*********************************************************************************************************************************** */

	template<class Type>
	class sorter {
	public:
		std::vector<int> perm;
	//	sortperm: sorts the perm vector in such a way that 
	//	-	list[perm[i]] >= list[perm[i']] for all i' >= i, if input sense is 'u'
	//	-	list[perm[i]] <= list[perm[i']] for all i' >= i, if input sense is 'l'
		void sortperm(std::vector<Type> list, char sense) {
			perm.resize(list.size());
			for (int i = 0; i < perm.size(); i++) {
				perm[i] = i;
			}
			if (sense == 'l') std::sort(perm.begin(), perm.end(), [&list](int a, int b) { return list[a] < list[b]; });
			if (sense == 'u') std::sort(perm.begin(), perm.end(), [&list](int a, int b) { return list[a] > list[b]; });
		}
	};
