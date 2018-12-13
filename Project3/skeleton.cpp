#include<limits>
#include <numeric>
#include<iostream>
#include "cfl_data.h"
#include "knapsack.h"
#include "time.h"

using namespace std;
/************************************************************************************************************************************
	Function prototype
*************************************************************************************************************************************/
// this line makes this function visible to the rest of the file.
bool isfeasible(vector<bool> &OpenPlant, vector<std::vector<int>> x, cfldata &cfl, cflsol &sol);
void heuristic(cfldata &cfl, vector<std::vector<int>> &x, vector<bool> &OpenPlant, vector<double> &residualCapacity, vector<bool> &Custom, vector<sorter<double>> sortCustCosts, cflsol &sol);;
void kprepair(cfldata &cfl, vector<std::vector<int>> &x, vector<bool> &OpenPlant, vector<double> &residualCapacity, vector<bool> &Custom, vector<sorter<double>> sortCustCosts);
void repaire(cfldata &cfl, vector<std::vector<int>> &x, vector<bool> &OpenPlant, vector<double> &residualCapacity, vector<bool> &Custom, vector<sorter<double>> sortCustCosts, cflsol &sol);
void gapheu(cfldata &cfl, std::vector<bool> &OpenPlant, cflsol &sol);
void Lagrangean_solver(cfldata &cfl, cflsol &sol_star);

#define INFINITY DBL_MAX

// NOTE: "Plant" and "Facility" are used as synonims below

// *****************************************************************************************************************************************
// FUNCTION SOLVE SUBPROBLEM: 
// *****************************************************************************************************************************************
//  * Called |M| times at each iteration to solve the Lagrangean subproblem IP_j(u) of each facility j
// *****************************************************************************************************************************************


bool isfeasible(vector<bool> &OpenPlant, vector<std::vector<int>> x, cfldata &cfl, cflsol &sol) {
	bool isfea = false;
	for (int i = 0; i < cfl.nCustom; i++) {
		int isassigned = 0;
		for (int j = 0; j < cfl.nPlants; j++) {
			if (x[i][j] = 1) {
				isassigned += 1;
				sol.CusPlant[i] = 1;
				sol.cost += cfl.CustCosts[i][j];
			}
		}
		if (isassigned == 1) {
			isfea = true;
		}
		else {
			return false;
		}
	}
	if (isfea) {
		for (int j = 0; j < cfl.nPlants; j++) {
			if (OpenPlant[j]) {
				sol.cost += cfl.OpenCosts[j];
			}
		}
	}
	else {
		// this solution is not feasible, reset the solution;
		cflsol solt(cfl);
		sol = solt;
	}
	return isfea;
}


// *****************************************************************************************************************************************
// FUNCTION HEURISTIC: 
// *****************************************************************************************************************************************
//  * Called after Lagrangean problem IP(u) has been solved
//	* Implements one of the heuristics described in the project instructions and lecture. Returns a feasible solution and 
//		its cost ub, if one is found
// *****************************************************************************************************************************************

void heuristic(cfldata &cfl, vector<std::vector<int>> &x, vector<bool> &OpenPlant, vector<double> &residualCapacity, vector<bool> &Custom, vector<sorter<double>> sortCustCosts, cflsol &sol){
	// see lecture slides 03_cfl.pdf
	// GAP heuristic or Lagrangean repair heuristic?
	// std::cout << "-----------------Start: Information printed from function heuristic()------------------" << std::endl;
	std::cout << "Do GAP heuristic:" << std::endl;
	//repaire(cfl, x, OpenPlant, residualCapacity, Custom, sortCustCosts, sol);
	gapheu(cfl, OpenPlant, sol);
	// std::cout << "-----------------End: Information printed from function heuristic()--------------------" << std::endl;
}

void repaire(cfldata &cfl, vector<std::vector<int>> &x, vector<bool> &OpenPlant, vector<double> &residualCapacity, vector<bool> &Custom, vector<sorter<double>> sortCustCosts, cflsol &sol){
	int n = 0;
	vector<int> v; // [2,3,8] customer i is assigned to facility 2, 3 and 8 
	int bestj = 0;
	for (int i = 0; i < cfl.nCustom; i++) {
		for (int j = 0; j < cfl.nPlants; j++) {
			if (x[i][j] == 1) {
				n++;
				v.push_back(j);
				if (n > 1) {
					if (cfl.CustCosts[i][bestj] > cfl.CustCosts[i][j]) {
						x[i][bestj] = 0;
						bestj = j;
					}
					else {
						x[i][j] = 0;
					}
				}
				else {
					bestj = j;
				}
			}
		}
		// TODO: problem in above code: because some customer has been removed from facility, some facility maybe closed,
		// need to updata OpenPlant;

		// customer i is not assigned to any facility
		kprepair(cfl, x, OpenPlant, residualCapacity, Custom, sortCustCosts);
	 }
}

void kprepair(cfldata &cfl, vector<std::vector<int>> &x, vector<bool> &OpenPlant, vector<double> &residualCapacity, vector<bool> &Custom, vector<sorter<double>> sortCustCosts) {

	// Knapsack heuristic - repair x -parameter: x, sortCustCosts, residualCapacity, OpenPlant, Custom
	sorter<double> sortOpenCosts;
	sortOpenCosts.sortperm(cfl.OpenCosts, 'l');
	for (int i = 0; i < cfl.nCustom; i++) {
		bool tem = false;
		if (!Custom[i])
			sortCustCosts[i].sortperm(cfl.CustCosts[i], 'l');		//Sort serve cost of customer i for increasing opening cost;
		for (int j = 0; j < cfl.nPlants; j++) {
			if (OpenPlant[sortCustCosts[i].perm[j]] && residualCapacity[sortCustCosts[i].perm[j]] >= cfl.CustDeman[i] && !tem) {
				x[i][sortCustCosts[i].perm[j]] = 1;
				residualCapacity[sortCustCosts[i].perm[j]] -= cfl.CustDeman[i];
				Custom[i] = true;
				tem = true;
				break;
			}
		}
		if (!tem) {
			for (int j = 0; j < cfl.nPlants; j++) {
				if (!OpenPlant[sortOpenCosts.perm[j]]) {
					OpenPlant[sortOpenCosts.perm[j]] = true;
					x[i][sortOpenCosts.perm[j]] = 1;
					residualCapacity[sortOpenCosts.perm[j]] -= cfl.CustDeman[i];
					Custom[i] = true;
					break;
				}
			}
		}
	}
}

// A function that executes the GAP heuristic using the subset of open facilities 
// found by cheapest_facility_subset().
void gapheu(cfldata &cfl, std::vector<bool> &OpenPlant, cflsol &sol) {
	
	bool feas = true;
	/*
	std::cout << "Check initial value of variables needed:" << std::endl; // z solution cost
	std::cout << sol.cost << std::endl; // z solution cost
	*/
	std::vector<double> Qbar(cfl.nPlants);
	for (int j = 0; j < cfl.nPlants; j++) {
		Qbar[j] = cfl.PlantCapa[j];  //Initial Qbar the rest of Capacity of facilities
	}

	std::vector<int> U(cfl.nCustom);
	for (int i = 0; i < cfl.nCustom; i++) {
		U[i] = i;
	}

	//std::cout << "-----------------Start: Information printed from function gapheu()------------------" << std::endl;
	while (!U.empty() && feas) {
		double delta_star = -std::numeric_limits<double>::infinity();
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
		if (i_star != -1 && j_star != -1) {
			sol.CusPlant[i_star] = j_star;
			sol.cost += cfl.CustCosts[i_star][j_star];
			Qbar[j_star] -= cfl.CustDeman[i_star];
			for (int i = 0; i < U.size(); i++) {
				if (U[i] == i_star) {
					U[i] = U[U.size() - 1];
					break;
				}
			}
			U.resize(U.size() - 1);
			//std::cout << i_star << " is assigned to " << j_star << " U size:" << U.size() << std::endl;
		}
	}
	//std::cout << "-----------------End: Information printed from function gapheu()--------------------" << std::endl;
}


// *****************************************************************************************************************************************
// LAGRANGEAN ALGORITHM FUNCTION: 
// *****************************************************************************************************************************************

void Lagrangean_solver(cfldata &cfl, cflsol &sol_star) {
	// Define the data structures needed. 
	int MAX_t =50;						// Should be [800,2000]
	int t = 0 ;
	std::vector<std::vector<double>> u(MAX_t+1, std::vector<double>(cfl.nCustom,0)); // the Lagrangean multipliers 
																				   // Initialize all the Lagrangean multipliers u_i equal to 0 
	std::vector<double> LB(MAX_t, 0);             // For remembering the low bound in each iteration of subgradient algo.
											   // this can computed by adding dpknaps.totalModifedcost and u[t].  
											   //(But we don't store dpknaps.totalModifedcost of each itreation, and we only need this value,
											   // don't need total solution.)
	std::vector<double> z_ub(MAX_t, 0);           // For remembering the upper(primal) bound in each iteration of subgradient algo. 
											   // But this is also stored in sol[t].cost.   SO maybe we don't need this??
	std::vector<cflsol> sol;

	// Initialize the best dual and primal bound to -Inf and +Inf, respectively
	double LB_star = -DBL_MAX;             // LB i is the best dual bound so far
	double z_ub_star = DBL_MAX;            // z_ub is the best know primal bound
	std::vector<double> u_star(cfl.nCustom);

	// Set the scaling factor alpha for the step size, e.g., alpha. 
	// If using the step update rule 1 then theta = lambda 
	double alpha = 0.996;            // update rule 1, should try [0.5,0.999]
	double theta = 2.0;			  // update rule 1, should try [1.0,2.0]
	double lambda = 2.0;
	double epsilon = 0.0002;

	// Compute a trivial upper bound ub: Sum of all opening costs plus the sum over i of the largest c_(ij) 
	std::vector<sorter<double>> sortCustCosts(cfl.nCustom);
	for (int i = 0; i < cfl.nCustom; i++) {
		sortCustCosts[i].sortperm(cfl.CustCosts[i], 'l');		//Sort serve cost of customer i for increasing opening cost;
		z_ub[t] = z_ub[t] + cfl.CustCosts[i][sortCustCosts[i].perm[cfl.nPlants-1]];
	}
	z_ub[t] = z_ub[t] + std::accumulate(cfl.OpenCosts.begin(), cfl.OpenCosts.end(), 0);
	if (z_ub[t] < z_ub_star) {
		z_ub_star = z_ub[t];
	}
	std::cout << "trivial z_ub_star:"<< z_ub_star << std::endl;
	
	//  MAIN LOOP OF SUBGRADIENT ALGO:
	sorter<double> sortPlantCapa;
	sortPlantCapa.sortperm(cfl.PlantCapa, 'u');
	dpknap kp(cfl.PlantCapa[sortPlantCapa.perm[0]], cfl.nCustom);
	while (t < MAX_t -1 || theta < epsilon)
	{
		// For each facility j: do kp; 
		// At the same time construct solution x[i][j], compute the 'OpenPlant', 'Custom', 'residualCapacity' and partial LB[t] (i.e. partial totalModifiedCost) 
		std::cout << "--------------Start: Information printed from iteration:" << t << "-----------------------" << std::endl;
		double totalModifiedCost = 0;
		std::vector<std::vector<int>> x(cfl.nCustom, std::vector<int>(cfl.nPlants, 0));
		std::vector<bool> OpenPlant(cfl.nPlants, false); // compute the OpenPlant
		std::vector<bool> Custom(cfl.nCustom, false);
		std::vector<double> residualCapacity(cfl.nPlants, 0);
		residualCapacity = cfl.PlantCapa;
		cflsol sol_t(cfl);
		sol.push_back(sol_t); // the i-th iteration solution
		
		for (int j = 0; j < cfl.nPlants; j++) {
			/* 
				Prepare the data for the j-th Lagrangean subproblem (knapsack problem): Set the profit of each item i
				as - (c_(ij) - u_i) and its weight as q(i). The Knapsack capacity is Q(j).
				!!!! NOTE: Only items i with *positive profit* and with q(i) <= Q(j) are cosidered when solving the knapsak problem. !!!!
				* Solve the j-th knapsack problem, get its optimal cost z(j). Construct the optimal solution x of corresponding subproblem.
					Remember to map each knapsack item into the corresponding customer:
					x_ij should be 1 if the item corresponding to customer i is in the optimal solution of the knapsack problem.
				* Store the subproblem solution and optimal cost of the subproblem
				kp.solve(Q, ncustomers, profits, demands);
			*/
			double pro;
			int nCustoKp = 0;
			std::vector<int> indexCustoKp;  // eg. custom 2,4,6 into kp 
			std::vector<double> CustDemantoKp; // 0, 1, 2 --- 2,4,6
			std::vector<double> profite; // 0, 1, 2 --- 2,4,6

			for (int k = 0; k < cfl.nCustom; k++) {
				pro = -(cfl.CustCosts[k][j] - u[t][k]);
				if (pro > 0) {
					nCustoKp++;
					indexCustoKp.push_back(k);
					profite.push_back(pro);
					CustDemantoKp.push_back(cfl.CustDeman[k]);
				}
			}
			if (nCustoKp != 0) {
				kp.solve(cfl.PlantCapa[j], nCustoKp, profite, CustDemantoKp);
				//solve_subproblem(kp, cfl.PlantCapa[j], cfl.nCustom, profite, cfl.CustDeman);
				if (kp.solutionOpt && cfl.OpenCosts[j] - kp.solutionCost < 0) {
					totalModifiedCost = totalModifiedCost + cfl.OpenCosts[j] - kp.solutionCost;
					//std::cout << "totalcost" << totalModifiedCost << std::endl;
					OpenPlant[j] = true;
					for (int i = 0; i < kp.n_in_solution; i++) {
						//std::cout <<" kp.solution[i] indexCustoKp[kp.solution[i]]"<< kp.solution[i] << indexCustoKp[kp.solution[i]] << std::endl;
						x[indexCustoKp[kp.solution[i]]][j] = 1;
						Custom[indexCustoKp[i]] = true;
						residualCapacity[j] -= cfl.CustDeman[indexCustoKp[kp.solution[i]]];
					}
				}
			}
			else {
				std::cout << "not do kp" << std::endl;
			}
		}
		// End Loop -- For each facility j: do kp

		// If the x is feasible, can return, it should also be the optimal.
		if (isfeasible(OpenPlant, x, cfl, sol[t])) {
			std::cout << "After kp, solution is feasible!";
			getchar();
			LB[t] = totalModifiedCost + std::accumulate(u[t].begin(), u[t].end(), 0);
			LB_star = LB[t];
			z_ub_star = sol[t].cost;
			break;
		}
		// Compute the subgradient s using the solution x(u). For each i, s_i is 1 minus the sum over j of x(u)_ij
		std::vector<double> s(cfl.nCustom);
		for (int i = 0; i < cfl.nCustom; i++) {
			s[i] = 1 - std::accumulate(x[i].begin(), x[i].end(), 0);
		}

		// Compute L(u) as the sum over j of min[0, f(j) - z(j)] plus the sum of all Lagrangean multipliers u_i. Here, f(j) denotes the opening cost of 
		// facility j and z(j) is the optimal cost of the j-th knapsack problem
		LB[t] = totalModifiedCost + std::accumulate(u[t].begin(), u[t].end(), 0);

		// Call function HEURISTIC to execute the Lagrangean Heuristic: If a feasible solution is found and its cost is better than the best primal bound 
		// ub, update the best primal bound ub and store the solution

		kprepair(cfl, x, OpenPlant, residualCapacity, Custom, sortCustCosts);
		heuristic(cfl, x, OpenPlant, residualCapacity, Custom, sortCustCosts, sol[t]);
		sol[t].testfeas(cfl);
		if (sol[t].feas && sol[t].cost < z_ub_star) {
			z_ub_star = sol[t].cost;
			sol_star = sol[t];
		}

		// If L(u) is greater than the best dual bound then
		//	* update the best dual bound
		//	* store the current Lagrangean multipliers u in the vector u* of best Lagrangean multipliers
		// Endif 
		if (LB[t] > LB_star) {
			LB_star = LB[t];
			u_star = u[t];
		}

		// If using step update rule 2 and L(u) has not improved for a prescribed number of consecutive iterations, update lambda = lambda * alpha 
		int number = 5;
		int temp = t;
		bool signal = true;
		if (t > 20) {
			while (number--) {
				if (LB[t] != LB[t - number]) {
					signal = false;
					break;
				}
			}
			if (signal) {
				lambda = lambda * alpha;
			}
		}
		// If using the update rule 2
		//	compute the squared norm S of the subgradient as the sum over i of s(i)*s(i)
		//  compute the step size theta = lambda * (ub - L(u)) / S;
		double S = 0;
		for (int i = 0; i < s.size(); i++) {
			S = S + s[i] * s[i];
		}
		//std::cout << "S:" << S << std::endl;
		theta = lambda * (z_ub_star - LB[t]) / S;
		
		// 	If using the update rule 1
		//   compute the step size theta = theta * alpha;
		// theta = theta * alpha;

		//  Update each Lagrangean multipliers: u(i) := u(i) + theta * s(i) for each customer i
		for (int i = 0; i < cfl.nCustom; i++) {
			u[t+1][i] = u[t][i] + theta * s[i];
			//std::cout << u[t+1][i] << " ******************* ";
		}
		std::cout << "--------------End: Information printed from iteration:" << t << "-----------------------" << std::endl;
		t = t + 1;

		std::cout << "theta:" << theta << std::endl;
		std::cout << "z_ub_star: " << z_ub_star << " LB_star: " << LB_star  << std::endl;

		// 	Terminate the loop if 
		//	* a max number of iterations have been done
		//	* theta is smaller than an epsilon (e.g., 0.0002)
		//  * the opimality conditions for x(u) are verified    ????
		//  END MAIN LOOP OF SUBGRADIENT ALGO
	}
	
	//  Return the best primal and dual bound, the best feasible solution, the vector u*, the total execution time
	std::cout << "-------------------Start: Information printed from main function-----------------------------" << std::endl;
	std::cout << "z_ub_star: " << z_ub_star << " LB_star: " << LB_star << std::endl;
	std::cout << "" << std::endl;
	std::cout << "-------------------------------------------" << std::endl;


	std::cout << "sol_star:" << " ";
	sol_star.print(cfl, "filename");
	std::cout << "" << std::endl;
	std::cout << "-------------------------------------------" << std::endl;

	std::cout << "LB:" <<std::endl;
	for (int i = 0; i < MAX_t; i++) {
		std::cout << LB[i] << " ";
	}
	std::cout << "" << std::endl;
	std::cout << "-------------------------------------------" << std::endl;

	std::cout << "u_star:" << " ";
	for (int i = 0; i < cfl.nCustom; i++) {
		std::cout << u_star[i] << " ";
	}
	std::cout << "" << std::endl;
	std::cout << "-------------------End: Information printed  from main function------------------------------" << std::endl;
}


int main() {
	double time = 0;
	clock_t start, end;
	start = clock();
	cfldata cfl;
	bool quit = !(cfl.read("p2000-2000-89b.dat"));  // TODO: Change to read parameter from system args[1]
	std::cout << quit << std::endl;
	if (quit) return 0;

	cflsol sol(cfl);

	std::cout << "-------------------Start: Information printed from main function-----------------------------" << std::endl;
	std::cout << "Capacitated Facility Location Problem Name:" << cfl.probname << std::endl;
	std::cout << "The number of facilities:" << cfl.nPlants << std::endl;
	std::cout << "The number of customers:" << cfl.nCustom << std::endl;
	std::cout << "-------------------End: Information printed  from main function------------------------------" << std::endl;

	// Call the Lagrangean_solver()
	Lagrangean_solver(cfl, sol);
	end = clock();
	time = (end - start) / CLOCKS_PER_SEC;
	cout << " time: " << time << endl;
	getchar();
	getchar();
}




