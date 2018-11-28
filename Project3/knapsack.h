#include <vector>

class dpknap {
public:
	double solutionCost;                        // The solution cost
	std::vector<int> solution;                  // List of customers in the optimal solution
	int n_in_solution;                          // Number of customers in the optimal solution
	bool solutionOpt;
	std::vector<std::vector<double>> F;         // The DP table
	long maxcapacity;
	int maxitems;

	//  Constructor of a dpknap object. You must give in input:
	//	- the maximum capacity of any instance that will be solved (maxcapacity)
	//  - the maximum # of customers of any instance that will be solved

	dpknap(int maxcapacity, int maxitems) : maxcapacity(maxcapacity), maxitems(maxitems) {

		F.resize(maxcapacity + 1, std::vector<double>(maxitems+1));
		solution.resize(maxitems+1);
		n_in_solution = 0;
		solutionCost = 0;
		solutionOpt = false;
	}

	//  Function to solve a knapsack instance 
	
	//  INPUT 
	//	- the knapsack capacity (Q)
	//  - the # of customers (ncustomers)
	//  - the vector of customer profits (profits)
	//  - the vector of customer demands (demands)
	//    NOTE: The number of items and knapsack capacity passed in input cannot exceed the values 'maxitems' 
	//    and 'maxcapacity', respectively, passed in input to the object constructor
	
	//  OUTPUT
	//	At termination, the solution and its cost are obtained as follows	
	//		- solutionCost stores the cost of best solution found
	//		- solutionOpt is true if the solution found is proved optimal
	//		- n_in_solution stores the number of customers in the solution found
	//		- solution[0,...,n_in_solution] stores the subset of customers in the solution found
	
	//	 NOTE1: Returns TRUE wihout computing any solution if the # of items or capacity in input exceeds  
	//	 the maximum passed in input to the object constructor
	
	//   NOTE2: It is assumed that all input customer profits and weigths are positive 

	bool solve(int Q, int ncustomers, std::vector<double> profits, std::vector<double> demands) {

		solutionCost = 0;
		solutionOpt = false;
		n_in_solution = 0;

		if (ncustomers > maxitems || Q > maxcapacity) {
			return true;
		}

		int ibeg = 0;
		while (demands[ibeg] > Q) {
			ibeg++;
		}

		if (ibeg > ncustomers - 1) {
			return false;
		}
		for (int q = 0; q <demands[ibeg]; q++) {
			F[q][ibeg] = 0;
		}
		for (int q = demands[ibeg]; q <= Q; q++) {
			F[q][ibeg] = profits[ibeg];
		}
		for (int i = ibeg+1; i < ncustomers; i++) {
			for (int q = 0; q <= Q; q++) {
				F[q][i] = F[q][i - 1];
				if (q >= demands[i]) {
					if (F[q - demands[i]][i - 1] + profits[i] > F[q][i]) {
						F[q][i] = F[q - demands[i]][i - 1] + profits[i];
					}
				}
			}
		}
		solutionCost = F[Q][ncustomers - 1];
		solutionOpt = true;

		int load = Q;
		double profit = 0;
		int i = ncustomers - 1;
		while (i > ibeg) {
			if (F[load][i] > F[load][i - 1]) {
				solution[n_in_solution++] = i;
				load -= demands[i];
				profit += profits[i];
			}
			i--;
		}
		if (F[load][i] > 0) {
			solution[n_in_solution++] = i;
			load -= demands[i];
			profit += profits[i];
		}
		return false;
	}
};

// =================================================================================================================================================================
// EXAMPLE USAGE:
// =================================================================================================================================================================
// Call the constructor to create a knapsack solver object (kp). 
// The input parameters are the maximum knpasack capacity (Q_max) and the maximum number of items (n_max)
// For the project, Q_max is the largest facility capacity, and n_max is the number of customers

// dpknap kp(Q_max, n_max);

// Call the knpasack solve function to solve a knapsack instance.
// The input parameters are the capacity of the knpasack (Q), the number of items (n), the vector of item profits (profits), and the vector of item weights (weights)

// kp.solve(Q, n, profits, weights);

// If all goes well, after the call to kp.solve() the knapsack solution will be stored in the internal vector kp.solution
// The number of items in the optimal solution is stored in the internal variable kp.n_in_solution

