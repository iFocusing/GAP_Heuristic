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
bool isfeasible(vector<bool> &OpenPlant, vector<vector<int>> x, cfldata &cfl, cflsol &sol);
void heuristic(cfldata &cfl, vector<vector<int>> &x, vector<bool> &OpenPlant, vector<double> &residualCapacity, vector<bool> &Custom, vector<sorter<double>> sortCustCosts, cflsol &sol);;
void repaire(cfldata &cfl, vector<vector<int>> &x, vector<bool> &OpenPlant, vector<double> &residualCapacity, vector<bool> &Custom, vector<sorter<double>> sortCustCosts);
void gapheu(cfldata &cfl, vector<bool> &OpenPlant, cflsol &sol);
void Lagrangean_solver(cfldata &cfl, cflsol &sol_star, int MAX_t, int number, int updaterule);
void gapheu_qq(cfldata &cfl, vector<int> &OpenPlant, cflsol &sol);

#define INFINITY DBL_MAX

// Folowing isfeasible function is to check if the solution after doing sub_problem
// (i.e. after doing kp for all facilities)
bool isfeasible(vector<bool> &OpenPlant, vector<vector<int>> x, cfldata &cfl, cflsol &sol) {
	bool isfea = false;
	for (int i = 0; i < cfl.nCustom; i++) {
		int isassigned = 0;
		for (int j = 0; j < cfl.nPlants; j++) {
			if (x[i][j] = 1) {
				isassigned += 1;
				sol.CusPlant[i] = j;
			}
		}
		if (isassigned == 1) {
			isfea = true;
		}
		else {
			return false;
		}
	}
	return sol.testfeas(cfl);
}

// *****************************************************************************************************************************************
// FUNCTION HEURISTIC: 
// *****************************************************************************************************************************************
//  * Called after Lagrangean problem IP(u) has been solved
//	* Implements one of the heuristics described in the project instructions and lecture. Returns a feasible solution and 
//		its cost ub, if one is found
// *****************************************************************************************************************************************

void heuristic(cfldata &cfl, vector<vector<int>> &x, vector<bool> &OpenPlant, vector<double> &residualCapacity, vector<bool> &Custom, vector<sorter<double>> sortCustCosts, cflsol &sol){
	// see lecture slides 03_cfl.pdf
	// GAP heuristic or Lagrangean repair heuristic
	repaire(cfl, x, OpenPlant, residualCapacity, Custom, sortCustCosts);
	gapheu(cfl, OpenPlant, sol);
	
}

void repaire(cfldata &cfl, vector<vector<int>> &x, vector<bool> &OpenPlant, vector<double> &residualCapacity, vector<bool> &Custom, vector<sorter<double>> sortCustCosts){
	cout << "Do repaire heuristic:" << endl;
	// following is step 1,2 in the lecture;
	for (int i = 0; i < cfl.nCustom; i++) {
		int n = 0;
		int bestj;
		for (int j = 0; j < cfl.nPlants; j++) {
			if (x[i][j] == 1) {
				n++;
				if (n > 1) {
					if (cfl.CustCosts[i][bestj] > cfl.CustCosts[i][j]) {
						x[i][bestj] = 0;
						residualCapacity[bestj] -= cfl.CustDeman[i];
						bestj = j;
						x[i][bestj] = 1;
					}
					else {
						x[i][j] = 0;
						residualCapacity[j] -= cfl.CustDeman[i];
					}
				}
				else {
					bestj = j;
				}
			}
		}
	}
	// modifyOpenPlant after doing step1, 2;
	for (int j = 0; j < cfl.nPlants; j++) {
		bool openOrNot = false;
		if (OpenPlant[j]) {
			for (int i = 0; i < cfl.nCustom; i++) {
				if (x[i][j] == 1) {
					openOrNot = true;
				}
			}
		}
		OpenPlant[j] = openOrNot;
	}
	// function kprepair does step 3 in the lecture;
	sorter<double> sortOpenCosts;
	sortOpenCosts.sortperm(cfl.OpenCosts, 'l');
	for (int i = 0; i < cfl.nCustom; i++) {
		bool tem = false;
		if (!Custom[i])
			sortCustCosts[i].sortperm(cfl.CustCosts[i], 'l');//Sort serve cost of customer i for increasing opening cost;
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
// found by repaire heristic or else.


void gapheu_qq(cfldata &cfl, vector<int> &OpenPlant, cflsol &sol) {
	vector<double> resiCap(cfl.nPlants, 0);
	vector<vector<int>> sortF(cfl.nCustom, vector<int>(cfl.nPlants));

	vector<int> cusPlant(cfl.nCustom, -1);
	vector<double> U(cfl.nCustom);
	bool feas = true;
	vector<int> F;

	for (int j = 0; j < cfl.nPlants; j++) {
		if (OpenPlant[j] == 1) {
			resiCap[j] = cfl.PlantCapa[j];
		}
	}

	sorter<double> sort;
	for (int i = 0; i < cfl.nCustom; i++) {
		U[i] = i;
		sort.sortperm(cfl.CustCosts[i], 'l');
		sortF[i] = sort.perm;
	}

	while (!U.empty() && feas) {
		// step2.a
		double delta_star = DBL_MIN;
		double j_star = -1;
		double i_star = -1;
		for (int i = 0; i < U.size(); i++) {
			vector<int>().swap(F);
			int j1 = -1;
			int j2 = -1;
			double delta = DBL_MAX;
			for (int j = 0; j < cfl.nPlants; j++) {
				if (resiCap[sortF[U[i]][j]] >= cfl.CustDeman[U[i]]) {
					F.push_back(sortF[U[i]][j]);
				}
			}
			if (F.size() == 0) {
				sol.feas = false;
				return;
			}
			else if (F.size() == 1) {
				j1 = F[0];
			}
			else {
				j1 = F[0];
				j2 = F[1];

				delta = cfl.CustCosts[U[i]][j2] - cfl.CustCosts[U[i]][j1];
			}
			if (delta > delta_star) {
				delta_star = delta;
				j_star = j1;
				i_star = i;
			}
		}
		cusPlant[U[i_star]] = j_star;
		resiCap[j_star] -= cfl.CustDeman[U[i_star]];
		vector<double>::iterator it = U.begin() + i_star;
		U.erase(it);

	}
	sol.CusPlant = cusPlant;
	sol.testfeas(cfl);
}



void gapheu(cfldata &cfl, vector<bool> &OpenPlant, cflsol &sol) {
	cout << "Do GAP heuristic:" << endl;
	bool feas = true;
	vector<double> Qbar(cfl.nPlants);
	for (int j = 0; j < cfl.nPlants; j++) {
		Qbar[j] = cfl.PlantCapa[j];  //Initial Qbar the rest of Capacity of facilities
	}
	vector<int> U(cfl.nCustom);
	for (int i = 0; i < cfl.nCustom; i++) {
		U[i] = i;
	}
	while (!U.empty() && feas) {
		double delta_star = -numeric_limits<double>::infinity();
		int j_star = -1;
		int i_star = -1; //we will assign customer i_star to facility j_star at the end of this while_loop 
		vector<int> F; //the set of index facilities that can serve customer i, temporly not stored
		for (int i = 0; i < U.size(); i++) {
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
				delta = numeric_limits<double>::infinity();
			}
			else {
				vector<double> CustCosts_F_i;
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
			F.clear();
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
		//	cout << i_star << " is assigned to " << j_star << " U size:" << U.size() << endl;
		}

	}
}


// *****************************************************************************************************************************************
// LAGRANGEAN ALGORITHM FUNCTION: 
// *****************************************************************************************************************************************

void Lagrangean_solver(cfldata &cfl, cflsol &sol_star, int MAX_t, int number, int updaterule) {
	// Define the data structures needed. 
	int t = 0 ;
	vector<vector<double>> u(MAX_t+1, vector<double>(cfl.nCustom,0)); // the Lagrangean multipliers 
																				   // Initialize all the Lagrangean multipliers u_i equal to 0 
	vector<double> LB(MAX_t, 0);             // For remembering the low bound in each iteration of subgradient algo.
											   // this can computed by adding dpknaps.totalModifedcost and u[t].  
											   //(But we don't store dpknaps.totalModifedcost of each itreation, and we only need this value,
											   // don't need total solution.)
	vector<double> z_ub(MAX_t, 0);           // For remembering the upper(primal) bound in each iteration of subgradient algo. 
											   // But this is also stored in sol[t].cost.   SO maybe we don't need this??
	vector<cflsol> sol;

	// Initialize the best dual and primal bound to -Inf and +Inf, respectively
	double LB_star = -DBL_MAX;             // LB i is the best dual bound so far
	double z_ub_star = DBL_MAX;            // z_ub is the best know primal bound
	vector<double> u_star(cfl.nCustom);

	// Set the scaling factor alpha for the step size, e.g., alpha. 
	// If using the step update rule 1 then theta = lambda 
	double alpha = 0.999;            // update rule 1, should try [0.5,0.999]
	double theta = 2.0;			  // update rule 1, should try [1.0,2.0]
	double lambda = 2.0;
	double epsilon = 0.0000002;

	// Compute a trivial upper bound ub: Sum of all opening costs plus the sum over i of the largest c_(ij) 
	vector<sorter<double>> sortCustCosts(cfl.nCustom);
	for (int i = 0; i < cfl.nCustom; i++) {
		sortCustCosts[i].sortperm(cfl.CustCosts[i], 'l');		//Sort serve cost of customer i for increasing opening cost;
		z_ub[t] = z_ub[t] + cfl.CustCosts[i][sortCustCosts[i].perm[cfl.nPlants-1]];
	}
	z_ub[t] = z_ub[t] + accumulate(cfl.OpenCosts.begin(), cfl.OpenCosts.end(), 0);
	if (z_ub[t] < z_ub_star) {
		z_ub_star = z_ub[t];
	}
	cout << "trivial z_ub_star:"<< z_ub_star << endl;
	
	//  MAIN LOOP OF SUBGRADIENT ALGO:
	sorter<double> sortPlantCapa;
	sortPlantCapa.sortperm(cfl.PlantCapa, 'u');
	dpknap kp(cfl.PlantCapa[sortPlantCapa.perm[0]], cfl.nCustom);
	vector<vector<int>> x(cfl.nCustom, vector<int>(cfl.nPlants, 0));
	vector<bool> OpenPlant(cfl.nPlants, false); // compute the OpenPlant
	vector<bool> Custom(cfl.nCustom, false);
	vector<double> residualCapacity(cfl.nPlants, 0);

	vector<int> indexCustoKp;  // eg. custom 2,4,6 into kp 
	vector<double> CustDemantoKp; // 0, 1, 2 --- 2,4,6
	vector<double> profite; // 0, 1, 2 --- 2,4,6

	while (t < MAX_t && theta > epsilon)
	{
		// For each facility j: do kp; 
		// At the same time construct solution x[i][j], compute the 'OpenPlant', 'Custom', 'residualCapacity' and partial LB[t] (i.e. partial totalModifiedCost) 
		cout << "--------------Start: Information printed from iteration:" << t << "-----------------------" << endl;
		// initial x
		for (int i = 0; i < cfl.nCustom; i++) {
			for (int j = 0; j < cfl.nPlants; j++) {
				x[i][j] = 0;
			}
		}
		// initial OpenPlant
		for (int j = 0; j < cfl.nPlants; j++) {
				OpenPlant[j] = false;
		}
		// initial Custom
		for (int i = 0; i < cfl.nCustom; i++) {
			Custom[i] = false;
		}
		// initial residualCapacity
		for (int j = 0; j < cfl.nPlants; j++) {
			residualCapacity[j] = cfl.PlantCapa[j];
		}
		double totalModifiedCost = 0;
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
				if (kp.solutionOpt && cfl.OpenCosts[j] - kp.solutionCost < 0) {
					totalModifiedCost = totalModifiedCost + cfl.OpenCosts[j] - kp.solutionCost;
					OpenPlant[j] = true;
					for (int i = 0; i < kp.n_in_solution; i++) {
						//cout <<" kp.solution[i] indexCustoKp[kp.solution[i]]"<< kp.solution[i] << indexCustoKp[kp.solution[i]] << endl;
						x[indexCustoKp[kp.solution[i]]][j] = 1;
						Custom[indexCustoKp[i]] = true;
						residualCapacity[j] -= cfl.CustDeman[indexCustoKp[kp.solution[i]]];
					}
				}
			}
			else {
				//cout << "Not do kp for this facility " << j <<", all profites are non-positive!" << endl;
			}

			indexCustoKp.clear();
			CustDemantoKp.clear();
			profite.clear();
		}
		// End Loop -- For each facility j: do kp

		// If the x is feasible, can return, it should also be the optimal.
		if (isfeasible(OpenPlant, x, cfl, sol[t])) {
			LB[t] = totalModifiedCost + accumulate(u[t].begin(), u[t].end(), 0);
			LB_star = LB[t];
			z_ub_star = sol[t].cost;
			break;
		}
		// Compute the subgradient s using the solution x(u). For each i, s_i is 1 minus the sum over j of x(u)_ij
		vector<double> s(cfl.nCustom);
		for (int i = 0; i < cfl.nCustom; i++) {
			s[i] = 1 - accumulate(x[i].begin(), x[i].end(), 0);
		}
		
		// Compute L(u) as the sum over j of min[0, f(j) - z(j)] plus the sum of all Lagrangean multipliers u_i. Here, f(j) denotes the opening cost of 
		// facility j and z(j) is the optimal cost of the j-th knapsack problem
		LB[t] = totalModifiedCost + accumulate(u[t].begin(), u[t].end(), 0);

		// Call function HEURISTIC to execute the Lagrangean Heuristic: If a feasible solution is found and its cost is better than the best primal bound 
		// ub, update the best primal bound ub and store the solution

		//kprepair(cfl, x, OpenPlant, residualCapacity, Custom, sortCustCosts);
		heuristic(cfl, x, OpenPlant, residualCapacity, Custom, sortCustCosts, sol[t]);
		sol[t].testfeas(cfl);
		if (sol[t].feas && sol[t].cost < z_ub_star) {
			cout << "after herustic solution is feasible?" << sol[t].feas << endl;
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

		
		
		if (updaterule == 1) {
			// 	If using the update rule 1
			//   compute the step size theta = theta * alpha;
			theta = theta * alpha;
		}
		else {
			// If using step update rule 2 and L(u) has not improved for a prescribed number of 
			// consecutive iterations, update lambda = lambda * alpha 
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
			double S = 0.0;
			for (int i = 0; i < s.size(); i++) {
				S = S + s[i] * s[i];
			}
			theta = lambda * (z_ub_star - LB[t]) / S;
		}

		//  Update each Lagrangean multipliers: u(i) := u(i) + theta * s(i) for each customer i
		for (int i = 0; i < cfl.nCustom; i++) {
			u[t+1][i] = u[t][i] + theta * s[i];
			//cout << u[t+1][i] << " ******************* ";
		}
		cout << "--------------End: Information printed from iteration:" << t << "-----------------------" << endl;
		t = t + 1;

		cout << "theta:" << theta << endl;
		cout << "z_ub_star: " << z_ub_star << " LB_star: " << LB_star  << endl;

		// 	Terminate the loop if 
		//	* a max number of iterations have been done
		//	* theta is smaller than an epsilon (e.g., 0.0002)
		//  * the opimality conditions for x(u) are verified    ????
		//  END MAIN LOOP OF SUBGRADIENT ALGO
	}
	
	//  Return the best primal and dual bound, the best feasible solution, the vector u*, the total execution time
	cout << "-------------------Start: Information printed from main function-----------------------------" << endl;
	cout << "z_ub_star: " << z_ub_star << " LB_star: " << LB_star << endl;
	cout << "" << endl;
	cout << "-------------------------------------------" << endl;


	cout << "sol_star:" << " ";
	sol_star.print(cfl, "filename");
	cout << "" << endl;
	cout << "-------------------------------------------" << endl;

	cout << "LB:" <<endl;
	for (int i = 0; i < MAX_t; i++) {
		cout << LB[i] << " ";
	}
	cout << "" << endl;
	cout << "-------------------------------------------" << endl;

	cout << "u_star:" << " ";
	for (int i = 0; i < cfl.nCustom; i++) {
		cout << u_star[i] << " ";
	}
	cout << "" << endl;
	cout << "-------------------End: Information printed  from main function------------------------------" << endl;
}


int main() {
	double time = 0;
	clock_t start, end;
	start = clock();
	cfldata cfl;
	int MAX_t = 500;	// Should be [800,2000]
	string filename;
	cout << "which file do you want to test:" << endl;
	cin >> filename;
	cout << "How many iterations do you want to test:" << endl;
	cin >> MAX_t;

	int number;
	int updaterule;
	cout << "Updatestep():" << endl;
	cin >> number;
	cout << "updaterule(1 or 2)" << endl;
	cin >> updaterule;


	bool quit = !(cfl.read(filename));  // TODO: Change to read parameter from system args[1]
	if (quit) return 0;
	

	cflsol sol(cfl);

	// cout << "-------------------Start: Information printed from main function-----------------------------" << endl;
	cout << "Capacitated Facility Location Problem Name:" << cfl.probname << endl;
	cout << "The number of facilities:" << cfl.nPlants << endl;
	cout << "The number of customers:" << cfl.nCustom << endl;
	//cout << "-------------------End: Information printed  from main function------------------------------" << endl;

	// Call the Lagrangean_solver()
	Lagrangean_solver(cfl, sol, MAX_t, number, updaterule);
	end = clock();
	time = (end - start) / CLOCKS_PER_SEC;
	cout << " time: " << time << endl;
	cout << "Click enter to quit!" << endl;
	getchar();
	getchar();
}




