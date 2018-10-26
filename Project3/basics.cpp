
//  ..................................................................................................................................
//  FILE INCLUSIONS
//  ..................................................................................................................................

	#include "cfl_data.h"		// our project data classes and general functions
	#include <limits>			// for int type limits
	#include <float.h>			// for float type limits




	void basics() {

//  ..................................................................................................................................
//  VARIABLE DECLARATIONS AND ASSIGNMENTS
//  ..................................................................................................................................

//  (1) declare some variables of different types
	int			integer;
	long		longint;
	long long	verylongint;
	short		shortint;
	double		fractional;
	char		character;

//  (2) assign to each declared variable the largest value it can store
	integer    = INT_MAX;
	longint    = LONG_MAX;
	verylongint = LLONG_MAX;
	shortint   = SHRT_MAX;
	fractional = DBL_MAX;

//  (2) initialize a char variable
	character = 'a';



//  ..................................................................................................................................
//  SOME LOGICAL OPERATORS AND OPERATIONS
//  ..................................................................................................................................

	bool a = true;
	bool b = false;

	bool a_AND_b = (a && b);						// a_AND_b is assigned value false
	bool a_OR_b  = (a || b);						// a_OR_b  is assigned value true
	bool NOT_b   = !b;								// NOT_b   is assigned value true

	int  n1 = 0;
	int  n2 = 1;

	bool n1_LARGER_than_n2 = (n1 > n2);				// n1_LARGER_than_n2 is assigned value false
	bool n1_EQUALS_n2 = (n1 == n2);					// n1_LARGER_than_n2 is assigned value false
	bool n1_LARGER_OR_EQUAL_than_n2 = (n1 <= n2);	// n1_LARGER_than_n2 is assigned value true
	bool n1_NOT_EQUAL_to_n2 = (n1 != n2);			// n1_LARGER_than_n2 is assigned value true



//  ..................................................................................................................................
//  SOME CONTROL LOOPS
//  ..................................................................................................................................

//  IF ELSE: executes the block of code inside the IF block if the bool variable condition is true, otherwise executes the ELSE block
	bool condition = (n1 > n2);
	if (condition) {
		// do something
	}
	else {
		// do something else
	}

//  FOR: executes the block of code inside the FOR LOOP a number nloop of times
	int nloop = 10;
	for (int i = 0; i < nloop; i++) {
		// do something
	}

//  WHILE: executes the block of code inside the WHILE LOOP as long as the bool variable loop is true
	bool loop = true;
	while (loop) {
		// do something
		// if (something) 
		// then loop = false	
	}



//  ..................................................................................................................................
//  SOME BASIC USAGE OF VECTORS
//  ..................................................................................................................................

//  (1) create a vector of integers of size 100 and initialize all its elements to be 0
	int size = 100;
	int value = 0;
	std::vector<int> test_vector(size, value);


//  (2) access the current size of the vector
	int actual_size = (int)test_vector.size();


//  (3) access element 10 of the vector, modify it to be 1
	int position = 10;
	int test_value = test_vector[position];	
	test_vector[position] = 1;

//  (5) create a 2d vector (matrix) of integers with n_rowws rows and n_cols columns
	int n_rows = 10;
	int n_cols = 5;
	std::vector<std::vector<int>> test_matrix(n_rows, std::vector<int>(n_cols));

//  (6) initialize all entries of the 2d vector to 0
	for (int i = 0; i < test_matrix.size(); i++) {
		for (int j = 0; j < test_matrix[i].size(); j++) {
			test_matrix[i][j] = 0;
		}
	}



//  ..................................................................................................................................
//  USAGE OF THE CUSTOM CLASS SORTER: TO OBTAIN A SORTED PERMUTATION OF THE ELEMENTS OF A VECTOR WITHOUT CHANGING THE VECTOR
//  NOTE: sorter IS A CUSTOM CLASS DEFIEND IN cfldata.h
//  ..................................................................................................................................

//  (1) create a sorter object 
	sorter<int> srt;

//	(2) using the sorter, get a permutation of test_vector elements arranged for non-increasing value
	srt.sortperm(test_vector, 'u');
//  srt.perm[i], i=0,1.., now stores the position of the i-th largest element of test_vector

//  (3) using the sorter, get a permutation of test_vector elements arranged for non-decreasing value
	srt.sortperm(test_vector, 'l');
//  srt.perm[i], i=0,1.., now stores the position of the i-th smallest element of test_vector




//  ..................................................................................................................................
//  SOME STRING OPERATIONS
//  ..................................................................................................................................

//  (1) create two strings test_string_1 and test_string_2 
	std::string test_string_1;
	test_string_1 = "test";

	std::string test_string_2;
	test_string_2 = "string";

//	(2) concatenate them into a third sting test_string_3
	std::string test_string_3;
	test_string_3 = test_string_1 + test_string_2;



//  ..................................................................................................................................
//  OUTPUT TO CONSOLE
//  ..................................................................................................................................

//  (1) output test_string_3 to console and move the cursor to the next line
	std::cout << "test_string_3 =" << test_string_3;
	std::cout << std::endl;

//  (2) output values of test_vector to console and move the cursor to the next line
	for (int i = 0; i < test_vector.size(); i++) {
		std::cout << "test_vector (" << i << ") =" << test_vector[i];
		std::cout << std::endl;
	}



//  ..................................................................................................................................
//  OPEN/READ_FROM/CLOSE INPUT FILE
//  ..................................................................................................................................

//  (1) open/close a file input stream. The file to open is assumed to be located in the same directory of the project file *.vcxproj

//  Create an INPUT stream attached to an input file named "file.dat"
	std::ifstream input("file.dat");
//  Check if the file was opened correctly
	if (input.is_open()) {
//		starts with current position at the beginning of the file
//		extracts from the file the data from the current position until the first white space and stores it in the string variable name
//		Moves current position right after the white space
		std::string name;
		input >> name;
//		extracts from the file the data from the current position until the first white space and stores it in the int variable number
//		Read a string from the file
		int number;
		input >> number;
//		...
//		close the input file when done
		input.close();
	}



//  ..................................................................................................................................
//  OPEN/WRITE_TO/CLOSE OUTPUT FILE
//  ..................................................................................................................................

//  (1) open/close a file output stream. The file is created in the same directory of the project file *.vcxproj

//  Create an OUTPUT stream attached to an output file named "file.dat"
	std::ofstream sol("file.dat");
//  Check if the file was opened correctly
	if (sol.is_open()) {
//		starts with current position at the beginning of the file
//		write to file, at current position, the string in quotes "CFL instance " (including spaces) 
		sol << "Cost "; 
//		move current position to the beginning of the next line
		sol << std::endl;
//		write to file, at current position, the data stored in variable cost
		double cost = 100.1;
		sol << cost;
//		...
//		close the input file when done
		input.close();
	}

	}