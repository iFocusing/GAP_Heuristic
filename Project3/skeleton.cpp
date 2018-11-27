
#include "cfl_data.h"
#include "knapsack.h"

#define INFINITY DBL_MAX

// NOTE: "Plant" and "Facility" are used as synonims below

// *****************************************************************************************************************************************
// FUNCTION SOLVE SUBPROBLEM: 
// *****************************************************************************************************************************************
//  * Called |M| times at each iteration to solve the Lagrangean subproblem IP_j(u) of each facility j
// *****************************************************************************************************************************************

void solve_subproblem( /* ... */ ) {

	//	* Prepare the data for the j-th Lagrangean subproblem (knapsack problem): Set the profit of each item i 
	//		as - (c_(ij) - u_i) and its weight as q(i). The Knapsack capacity is Q(j). 

	//	!!!! NOTE: Only items i with *positive profit* and with q(i) <= Q(j) are cosidered when solving the knapsak problem. !!!!
	
	//	* Solve the j-th knapsack problem, get its optimal cost z(j). Construct the optimal solution x of corresponding subproblem. 
	//		Remember to map each knapsack item into the corresponding customer: 
	//		x_ij should be 1 if the item corresponding to customer i is in the optimal solution of the knapsack problem. 
	
	//  * Store the subproblem solution and optimal cost of the subproblem

}


// *****************************************************************************************************************************************
// FUNCTION HEURISTIC: 
// *****************************************************************************************************************************************
//  * Called after Lagrangean problem IP(u) has been solved
//	* Implements one of the heuristics described in the project instructions and lecture. Returns a feasible solution and 
//		its cost ub, if one is found
// *****************************************************************************************************************************************

void heuristic( /* ... */ ) {

	// see lecture slides 03_cfl.pdf

}

// *****************************************************************************************************************************************
// LAGRANGEAN ALGORITHM FUNCTION: 
// *****************************************************************************************************************************************

void Lagrangean_solver( /* ... */ ) {

	// Define the data structures needed. 

	// Initialize the best dual and primal bound to -Inf and +Inf, respectively

	// Set the scaling factor alpha for the step size, e.g., alpha. 
	// If using the step update rule 1 then theta = lambda 

	// Compute a trivial upper bound ub: Sum of all opening costs plus the sum over i of the largest c_(ij) 
	// Initialize all the Lagrangean multipliers u_i equal to 0

	//  MAIN LOOP OF SUBGRADIENT ALGO:

	// For each facility j:
	//	Call function SOLVE SUBPROBLEM
	//	End Loop

	// Construct the solution x(u) obtained by putting together all the optimal subproblem solutions x_j

	// Compute L(u) as the sum over j of min[0, f(j) - z(j)] plus the sum of all Lagrangean multipliers u_i. Here, f(j) denotes the opening cost of 
	// facility j and z(j) is the optimal cost of the j-th knapsack problem

	// Compute the subgradient s using the solution x(u). For each i, s_i is 1 minus the sum over j of x(u)_ij

	// Call function HEURISTIC to execute the Lagrangean Heuristic: If a feasible solution is found and its cost is better than the best primal bound 
	// ub, update the best primal bound ub and store the solution

	// If L(u) is greater than the best dual bound then
	//	* update the best dual bound
	//	* store the current Lagrangean multipliers u in the vector u* of best Lagrangean multipliers
	// Endif 

	// If using step update rule 2 and L(u) has not improved for a prescribed number of consecutive iterations, update lambda = lambda * alpha 

	// If using the update rule 2
	//	compute the squared norm S of the subgradient as the sum over i of s(i)*s(i)
	//  compute the step size theta = lambda * (ub - L(u)) / S;

	// 	If using the update rule 1
	//   compute the step size theta = theta * alpha;

	//  Update each Lagrangean multipliers: u(i) := u(i) + theta * s(i) for each customer i

	// 	Terminate the loop if 
	//	* a max number of iterations have been done
	//	* theta is smaller than an epsilon (e.g., 0.0002)
	//  * the otimality conditions for x(u) are verified

	//  END MAIN LOOP OF SUBGRADIENT ALGO

	//  Return the best primal and dual bound, the best feasible solution, the vector u*, the total execution time

}