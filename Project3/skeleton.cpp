#include<limits>
#include <numeric>
#include<iostream>
#include "cfl_data.h"
#include "knapsack.h"

using namespace std;
/************************************************************************************************************************************
	Function prototype
*************************************************************************************************************************************/
// this line makes this function visible to the rest of the file.
void cheapest_facility_subset(cfldata &cfl, std::vector<bool> &OpenPlant);
void gapheu(cfldata &cfl, std::vector<bool> &OpenPlant, cflsol &sol);
dpknap solve_subproblem(long Q_max, int n_max, int Q, int ncustomers, std::vector<double> profits, std::vector<double> demands);
void heuristic(cfldata &cfl, cflsol &sol);
void gapheu(cfldata &cfl, std::vector<bool> &OpenPlant, cflsol &sol);
void Lagrangean_solver(cfldata &cfl, cflsol &sol_star);


#define INFINITY DBL_MAX

// NOTE: "Plant" and "Facility" are used as synonims below

// *****************************************************************************************************************************************
// FUNCTION SOLVE SUBPROBLEM: 
// *****************************************************************************************************************************************
//  * Called |M| times at each iteration to solve the Lagrangean subproblem IP_j(u) of each facility j
// *****************************************************************************************************************************************

dpknap solve_subproblem(long Q_max, int n_max, int Q, int ncustomers, std::vector<double> profits, std::vector<double> demands) {

	//	* Prepare the data for the j-th Lagrangean subproblem (knapsack problem): Set the profit of each item i 
	//		as - (c_(ij) - u_i) and its weight as q(i). The Knapsack capacity is Q(j). 

	//	!!!! NOTE: Only items i with *positive profit* and with q(i) <= Q(j) are cosidered when solving the knapsak problem. !!!!
	
	//	* Solve the j-th knapsack problem, get its optimal cost z(j). Construct the optimal solution x of corresponding subproblem. 
	//		Remember to map each knapsack item into the corresponding customer: 
	//		x_ij should be 1 if the item corresponding to customer i is in the optimal solution of the knapsack problem. 
	
	//  * Store the subproblem solution and optimal cost of the subproblem
	dpknap kp(Q_max, n_max);
	kp.solve(Q, ncustomers, profits, demands);

	std::cout << "-----------------Start: Information printed from function solve_subproblem()------------------" << std::endl;
	std::cout << "The number of items in the optimal solution: " << kp.n_in_solution;
	for (int i = 0; i < kp.solution.size()-1; i++) {
		std::cout << "Print the solution:" << kp.solution[i];
		std::cout << std::endl;
	}
	std::cout << "-----------------End: Information printed from function solve_subproblem()--------------------" << std::endl;
	return kp;
}


// *****************************************************************************************************************************************
// FUNCTION HEURISTIC: 
// *****************************************************************************************************************************************
//  * Called after Lagrangean problem IP(u) has been solved
//	* Implements one of the heuristics described in the project instructions and lecture. Returns a feasible solution and 
//		its cost ub, if one is found
// *****************************************************************************************************************************************

void heuristic(cfldata &cfl, cflsol &sol) {
	// see lecture slides 03_cfl.pdf
	// GAP heuristic or Lagrangean repair heuristic?
	std::cout << "-----------------Start: Information printed from function heuristic()------------------" << std::endl;
	std::vector<bool> OpenPlant(cfl.nPlants, false);
	cheapest_facility_subset(cfl, OpenPlant);
	std::cout << "Do GAP heuristic:" << std::endl;
	gapheu(cfl, OpenPlant, sol);
	std::cout << "-----------------End: Information printed from function heuristic()--------------------" << std::endl;
}


// A function that executes the GAP heuristic using the subset of open facilities 
// found by cheapest_facility_subset().
void gapheu(cfldata &cfl, std::vector<bool> &OpenPlant, cflsol &sol) {
	
	bool feas = true;
	std::cout << "Check initial value of variables needed:" << std::endl; // z solution cost
	std::cout << sol.cost << std::endl; // z solution cost
	std::vector<double> Qbar(cfl.nPlants);
	for (int j = 0; j < cfl.nPlants; j++) {
		Qbar[j] = cfl.PlantCapa[j];  //Initial Qbar the rest of Capacity of facilities
		std::cout << Qbar[j] << std::endl;
	}

	std::vector<int> U(cfl.nCustom);
	for (int i = 0; i < cfl.nCustom; i++) {
		U[i] = i;
	}

	std::cout << "-----------------Start: Information printed from function gapheu()------------------" << std::endl;
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
		std::cout << i_star << " is assigned to " << j_star << " U size:" << U.size() << std::endl;
	}
	std::cout << "-----------------End: Information printed from function gapheu()--------------------" << std::endl;
}


// A function that finds the cheapest subset of facilities that can supply greater or equal
void cheapest_facility_subset(cfldata &cfl, std::vector<bool> &OpenPlant) {
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

	for (int i = 0; i < cfl.nCustom; i++) {
		sumdemands += cfl.CustDeman[i];
	}

	for (int j = 0; j < cfl.nPlants; j++) {
		sumcapacity += cfl.PlantCapa[sort.perm[j]];
		OpenPlant[sort.perm[j]] = true;
		if (sumcapacity > sumdemands) {
			break;
		}
	}

	std::cout << "-------------------Start: Information printed from cheapest_facility_subset function-----------------------------" << std::endl;
	for (int i = 0; i < cfl.nPlants; i++) {
		std::cout << "If open? facility:" << i << ":" << OpenPlant[i] << " ";
	}
	std::cout << "-------------------End: Information printed  from cheapest_facility_subset function------------------------------" << std::endl;
}


// *****************************************************************************************************************************************
// LAGRANGEAN ALGORITHM FUNCTION: 
// *****************************************************************************************************************************************

void Lagrangean_solver(cfldata &cfl, cflsol &sol_star) {

	// Define the data structures needed. 
	double time = 0;					// Using the matirial provided by lecturer.
	int MAX_t = 8;						// Should be [800,2000]
	int t = 0 ;
	std::vector<std::vector<double>> u(MAX_t, std::vector<double>(cfl.nCustom,0)); // the Lagrangean multipliers 
																				   // Initialize all the Lagrangean multipliers u_i equal to 0 
	std::vector<double> LB(MAX_t);             // For remembering the low bound in each iteration of subgradient algo.
											   // this can computed by adding dpknaps.totalModifedcost and u[t].  
											   //(But we don't store dpknaps.totalModifedcost of each itreation, and we only need this value,
											   // don't need total solution.)
	std::vector<double> z_ub(MAX_t);           // For remembering the upper(primal) bound in each iteration of subgradient algo. 
											   // But this is also stored in sol[t].cost.   SO maybe we don't need this??

	std::vector<cflsol> sol;
	cflsol sol_t(cfl);
	sol.push_back(sol_t);
	
	// Initialize the best dual and primal bound to -Inf and +Inf, respectively
	double LB_star = -DBL_MAX;             // LB i is the best dual bound so far
	double z_ub_star = DBL_MAX;            // z_ub is the best know primal bound
	std::vector<double> u_star(cfl.nCustom);

	// Set the scaling factor alpha for the step size, e.g., alpha. 
	// If using the step update rule 1 then theta = lambda 
	float alpha_1 = 0.996;            // update rule 1, should try [0.5,0.999]
	float theta_1 = 2.0;			  // update rule 1, should try [1.0,2.0]
	float epsilon = 0.0002;

	// Compute a trivial upper bound ub: Sum of all opening costs plus the sum over i of the largest c_(ij) 
	sorter<double> sort;
	for (int i = 0; i < cfl.nCustom; i++) {
		sort.sortperm(cfl.CustCosts[i], 'u');		//Sort serve cost of customer i for decreasing opening cost;
		z_ub_star = z_ub_star + sort.perm[0];
	}
	z_ub_star = z_ub_star + std::accumulate(cfl.OpenCosts.begin(), cfl.OpenCosts.end(), 0);

	
	//  MAIN LOOP OF SUBGRADIENT ALGO:
	while (t <= MAX_t || theta_1 < epsilon)
	{
		// For each facility j:
		//	Call function SOLVE SUBPROBLEM
		//	End Loop

		// Construct the solution x(u) obtained by putting together all the optimal subproblem solutions x_j

		// Compute L(u) as the sum over j of min[0, f(j) - z(j)] plus the sum of all Lagrangean multipliers u_i. Here, f(j) denotes the opening cost of 
		// facility j and z(j) is the optimal cost of the j-th knapsack problem

		std::vector<double> profite(cfl.nCustom);
		double totalModifiedCost = 0;
		std::vector<std::vector<int>> x;
		for (int j = 0; j < cfl.nPlants; j++) {
			dpknap kp(cfl.PlantCapa[j], cfl.nCustom);
			for (int k = 0; k < cfl.nCustom; k++) {
				profite[k] = -(cfl.CustCosts[k][j] - u[t][k]);
			}
			kp = solve_subproblem(cfl.PlantCapa[j], cfl.nCustom, cfl.PlantCapa[j], cfl.nCustom, profite, cfl.CustDeman);
			if (cfl.OpenCosts[j] - kp.solutionCost < 0)
				totalModifiedCost = totalModifiedCost + cfl.OpenCosts[j] - kp.solutionCost;
			if (kp.solutionOpt) {
				for (int i = 0; i <= kp.n_in_solution; i++) {
					x[kp.solution[i]][j] = 1;
				}
			}
			else {
				std::cout << "--------------Start: Information printed from opereting an instance of class ipusol-----------------------" << std::endl;
				std::cout << "dp don't find a feasible solution for facility" << j << ", this means that the demands of \
								all customer is greater than the capacity of facility " << j << "!" << std::endl;
				std::cout << "--------------------End: Information printed  from main function--------------------------------------" << std::endl;
			}
		}

		
		LB[t] = totalModifiedCost + std::accumulate(u[t].begin(), u[t].end(), 0);

		// Compute the subgradient s using the solution x(u). For each i, s_i is 1 minus the sum over j of x(u)_ij
		std::vector<double> s(cfl.nCustom);
		for (int i = 0; i < cfl.nCustom; i++) {
			s[i] = 1 - std::accumulate(x[i].begin(), x[i].end(), 0);
		}

		// Call function HEURISTIC to execute the Lagrangean Heuristic: If a feasible solution is found and its cost is better than the best primal bound 
		// ub, update the best primal bound ub and store the solution
		heuristic(cfl, sol[t]);
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

		// If using the update rule 2
		//	compute the squared norm S of the subgradient as the sum over i of s(i)*s(i)
		//  compute the step size theta = lambda * (ub - L(u)) / S;

		// 	If using the update rule 1
		//   compute the step size theta = theta * alpha;
		theta_1 = theta_1 * alpha_1;

		//  Update each Lagrangean multipliers: u(i) := u(i) + theta * s(i) for each customer i
		for (int i = 0; i < cfl.nCustom; i++) {
			u[t][i] = u[t][i] + s[i];
		}
		t = t + 1;
		// 	Terminate the loop if 
		//	* a max number of iterations have been done
		//	* theta is smaller than an epsilon (e.g., 0.0002)
		//  * the opimality conditions for x(u) are verified    ????
		//  END MAIN LOOP OF SUBGRADIENT ALGO
	}

	 
	//  Return the best primal and dual bound, the best feasible solution, the vector u*, the total execution time
	std::cout << "-------------------Start: Information printed from main function-----------------------------" << std::endl;
	std::cout << z_ub_star << LB_star << "--sol_star--" << "---u_star---" << time << std::endl;
	std::cout << "-------------------End: Information printed  from main function------------------------------" << std::endl;
}


int main() {
	cfldata cfl;
	bool quit = !(cfl.read("cfl.dat"));  // TODO: Change to read parameter from system args[1]
	if (quit) return 0;

	cflsol sol(cfl);

	std::cout << "-------------------Start: Information printed from main function-----------------------------" << std::endl;
	std::cout << "Capacitated Facility Location Problem Name:" << cfl.probname << std::endl;
	std::cout << "The number of facilities:" << cfl.nPlants << std::endl;
	std::cout << "The number of customers:" << cfl.nCustom << std::endl;
	std::cout << "-------------------End: Information printed  from main function------------------------------" << std::endl;

	// Call the Lagrangean_solver()
	Lagrangean_solver(cfl, sol);
}

