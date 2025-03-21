// Huy Nguyen
// CS 6420
// Assignment #4
// Dr. Rague
// Due: 03/25/25
// Version: 1.0
// -----------------------------------------------------------------
// This program performs matrix multiplication using iterative,
// recursive, and Strassen's methods and compares their performance.
// -----------------------------------------------------------------
// Disclaimer:
// The requirement is to create function headers as followed
// double* iterative_mult(int n, double* A, double* B) => Matrix iterative_mul(Matrix &m)
// double* recursive_mult(int n, double* A, double* B) => Matrix recursive_mult(Matrix &m)
// double* strassen_mult(int n, double* A, double* B) => Matrix strassen_mul(Matrix &m)
// -----------------------------------------------------------------
// COMMENTS:
// **Reflect on the actual results you are getting and what you think is needed to make strassen_mult outperform the other two functions.**
// 1. So the expected result would be that the Strassen method would outperform the iterative and recursive methods.
// The order of runtime complexity for matrix multiplication is: Strassen < Recursive < Iterative
// That is because Strassen theoretically works under the assumption that only the multiplication is the most expensive operation. However, on a real machine, the intermediate addition, subtraction and multiplication recursive call require many intermediate sub-matrices creation, which is expensive; and the Strassen method is not cache-friendly, because it requires many sub-matrices to be created and destroyed in each recursive call. This makes the overhead of the Strassen method higher than the iterative method. 

// 2. Theoretically speaking, as the matrices get larger and larger, the overhead of the Strassen method will be less significant compared to the iterative method. The Strassen method will outperform the iterative method when the matrices are large enough. (Check strassen_over_iter.png), the chart depicts the ratio of Strassen over Iterative method. The ratio is still greater than 1, which means the Strassen method is still slower than the iterative method. However, the ratio is decreasing as the matrix size increases. This hints that the Strassen method will outperform the iterative method when the matrix size is large enough.
// Another way to improve the Strassen method is to use the Strassen method for large matrices and switch to the iterative method for small matrices. 
// (check comparison(strassen_iter_64).png), This chart depicts the runtime of Iterative, Recursive, and Strassen(Use Iterative for matrices leq 64) methods for square matrices size 2 to 2048. For matrices under 64, Strassen is the slow, but still in the same order of magnitude as the iterative method. For matrices larger than 64, the Strassen method starts to outperform the iterative method.

// runtime(org_strassen).png shows the runtime of iterative, recursive and original strassen.
// runtime(strassen_iter_64).png shows the runtime of iterative, recursive and strassen(use iterative for matrices leq 64)

// SUBMISSION:
// given_tests.txt contains the output of the test case for the given 2x2, 4x4, 10x10 matrices
// seven_tests.txt contains the output of the 7 different sizes of matrices 2, 4, 6, 8, 10, 13, 15

// Usage: g++ -std=c++11 matrix_mult_1D.cpp -o matrix_mult
// ./matrix_mult
// or
// make matrix_mult
// ./matrix_mult


// -----------------------------------------------------------------
// Compiler directives
#include <iostream>
#include <vector>
#include <iomanip>
#include <chrono>
// use pair
#include <utility>
#include <fstream>
#include <sstream>
#include <tuple>
#include <cmath>



// -----------------------------------------------------------------
// This class represents a matrix and provides methods for matrix
// operations including addition, subtraction, and multiplication.
//
// Version 1.0
// ----------------------------------------------------------
template <typename T> // use template for generic type T
class Matrix{
	private:
		T* matrix; // 1D array of type T(double for this assignment) to store the matrix column major order
		int rows;
		int cols;
	public:
		// Constructor to initialize the matrix with the given data and dimensions
		// arguments: rows, cols of the matrix, data to initialize the matrix
		Matrix(int rows, int cols, T* data){
			this->rows = rows;
			this->cols = cols;
			matrix = new T[rows*cols];
			for(int i=0; i<rows; i++){
				for(int j=0; j<cols; j++){
					// for row major order use : i*cols + j
					// for column major order use : i + j*rows
					matrix[i + j*rows] = data[i + j*cols];
				}
			}
		}

		// Constructor to initialize the matrix with zeros
		// arguments: rows, cols of the matrix

		Matrix(int rows, int cols){
			this->rows = rows;
			this->cols = cols;

			// Initialize matrix with zeros
			matrix = new T[rows*cols];

			for(int i=0; i<rows; i++){
				for(int j=0; j<cols; j++){
					matrix[i + j*rows] = 0;
				}
			}
		}
		
		// Get a submatrix of the current matrix
		// arguments: start_row, start_col, end_row, end_col of the submatrix
		// returns: submatrix of the current matrix
		Matrix get_sub_matrix(int start_row, int start_col, int end_row, int end_col){
			int sub_rows = end_row - start_row + 1;
			int sub_cols = end_col - start_col + 1;
			Matrix sub_matrix(sub_rows, sub_cols);

			for(int i=0; i<sub_rows; i++){
				for(int j=0; j<sub_cols; j++){
					sub_matrix(i, j) = (*this)(start_row + i, start_col + j);
				}
			}
			return sub_matrix;
		}
		//Getters
		int get_rows(){ // get the number of rows
			return this->rows;
		}

		int get_cols(){ // get the number of columns
			return this->cols;
		}

		//Setters
		void set_rows(int rows){ // set the number of rows
			this->rows = rows;
		}

		void set_cols(int cols){ // set the number of columns
			this->cols = cols;
		}

		// Get the matrix
		T*& get_matrix(){
			return matrix;
		}

		// Set operations A(i, j) value, where A is the matrix row i, column j
		T& operator()(int i, int j){
			// for row major order use : i*cols + j
			return matrix[i + j*rows];
		}

		// Matrix addition operator
		Matrix operator+(Matrix &m){
			if(this->rows != m.rows || this->cols != m.cols){
				std::cout << "Matrix addition not possible" << std::endl;
				return Matrix(0, 0);
			}
			Matrix result_matrix(this->rows, this->cols);

			for(int i=0; i<this->rows; i++){
				for(int j=0; j<this->cols; j++){
					//(*this) means get the object
					//(*this)(i,j) use the operator() to get the value at row i, column j
					result_matrix(i,j) = (*this)(i,j) + m(i,j);
				}
			}
			return result_matrix;
		}

		// Matrix subtraction
		Matrix operator-(Matrix &m){
			if(this->rows != m.rows || this->cols != m.cols){
				std::cout << "Matrix subtraction not possible" << std::endl;
				return Matrix(0, 0);
			}
			Matrix result_matrix(this->rows, this->cols);

			for(int i=0; i<this->rows; i++){
				for(int j=0; j<this->cols; j++){
					result_matrix(i, j) = (*this)(i,j) - m(i,j);
				}
			}
			return result_matrix;
		}



		//Print matrix in a readable format
		// arguments: none
		// returns: void
		// Data stored in column major order 1D array
		// matrix printed example:
		// | 1 2 3 |
		// | 4 5 6 |
		// | 7 8 9 |
		// stored in 1D array as [1, 4, 7, 2, 5, 8, 3, 6, 9]
		void print(){
			int longest_element = 0;

			// Find the longest element in the matrix
			// to format the matrix printing
			for(int i=0; i<this->rows; i++){
				for(int j=0; j<this->cols; j++){
					T current_element = (*this)(i,j);
						if(std::to_string((int)current_element).length() > longest_element){
						longest_element = std::to_string((int)current_element).length();
					}
				}
			}

			// Print the matrix
			for(int i=0; i<this->rows; i++){
				std::cout << "| ";
				for(int j=0; j<this->cols; j++){
					std::cout<< std::setw(longest_element) << (*this)(i,j) << " ";
				}
				std::cout << "|";
				std::cout << std::endl;
			}
		}

		//Print the matrices A, B and the result of their multiplication
		// arguments: Matrix A, Matrix B, Matrix result
		// returns: void
		// Example:
		// A
		// *
		// B
		// =
		// result
		void print_matrix_multiplication(Matrix &A, Matrix &B, Matrix &result){
			// This function prints the matrices A, B and the result of their multiplication
			// Assuming result matrix is already calculated from multiplying A and B
			A.print();
			std::cout  << "*" << std::endl;
			B.print();
			std::cout  << "=" << std::endl;
			result.print();

		}


		// Matrix multiplication
		// arguments: Matrix m, method to use for matrix multiplication, verbose to print the result
		// returns: time taken for matrix multiplication
		// The matrix object A is the left matrix and the matrix object m is the right matrix
		// example: A.multiply_matrix(B, "iter", true);
		// The above example multiplies matrix A with matrix B using the iterative method and prints the result
		T multiply_matrix(Matrix &m, std::string method="iter", bool verbose=true){
			if(method == "iter"){ // iterative method
				// Start the timer
				auto start = std::chrono::high_resolution_clock::now();

				// Multiply the matrices using the iterative method
				Matrix result = this->iterative_mul(m);

				// End the timer
				auto end = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed = end - start;

				// Print the time taken for matrix multiplication
				if (verbose){
					std::cout << "Time taken for matrix multiplication using iterative method: " << elapsed.count() << "s" << std::endl;
					this->print_matrix_multiplication(*this, m, result);
				}

				return elapsed.count();
			}
			else if(method == "recur"){
				// Start the timer
				auto start = std::chrono::high_resolution_clock::now();

				// Multiply the matrices using the recursive method
				Matrix result = this->recursive_mult(m);

				// End the timer
				auto end = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed = end - start;

				if (verbose) {
					std::cout << "Time taken for matrix multiplication using recursive method: " << elapsed.count() << "s" << std::endl;
					this->print_matrix_multiplication(*this, m, result);
				}
				return elapsed.count();
			}
			else if(method == "strassen"){
				// Start the timer
				auto start = std::chrono::high_resolution_clock::now();

				// Multiply the matrices using the Strassen method
				Matrix result = this->strassen_mul(m);

				// End the timer
				auto end = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed = end - start;

				if (verbose){
					std::cout << "Time taken for matrix multiplication using strassen method: " << elapsed.count() << "s" << std::endl;
					this->print_matrix_multiplication(*this, m, result);
				}
				return elapsed.count();
			}
			else{
				std::cout << "Invalid method" << std::endl;
			}
			return -1;
		}

		// The requirement is to create function headers as followed
		// double* iterative_mult(int n, double* A, double* B)

		// Matrix multiplication using the iterative method
		// This is called by the multiply_matrix method
		// but outside users can also call this method
		// arguments: Matrix m
		// returns: result matrix of the multiplication
		// The matrix object A is the left matrix and the matrix object m is the right matrix
		// example: A.iterative_mul(B);
		Matrix iterative_mul(Matrix &m){

			// Check if the matrices can be multiplied
			if(this->cols != m.rows){
				std::cout << "Matrix multiplication not possible" << std::endl;
				return Matrix(0, 0);
			}

			// Create a matrix to store the result of the multiplication
			Matrix result_matrix(this->rows, m.cols);

			// Multiply the matrices using the iterative method
			// Three nested loops to iterate over the rows and columns of the result matrix
			for(int i=0; i<this->rows; i++){
				for(int j=0; j<m.cols; j++){
					T sum = 0;
					for(int k=0; k<this->cols; k++){
						sum += (*this)(i, k) * m(k, j);
					}
					result_matrix(i,j) = sum;
				}
			}
			return result_matrix;
		}

		// Helper function for the recursive matrix multiplication
		// This function is called by the recursive_mult method
		private:
		void _recur_matrix(Matrix &A, Matrix &B, Matrix &C, int n, 
			std::pair<int, int> A_start={0,0}, std::pair<int, int> B_start={0,0}, std::pair<int, int> C_start={0,0}){

			if (n == 1){
				C(C_start.first, C_start.second) += A(A_start.first, A_start.second) * B(B_start.first, B_start.second);
			}
			else{
				// Divide the matrices into 4 submatrices
				std::pair<int, int> A11_start = A_start; // top left
				std::pair<int, int> A12_start = std::make_pair(A_start.first, A_start.second + n/2); // top right
				std::pair<int, int> A21_start = std::make_pair(A_start.first + n/2, A_start.second); // bottom left
				std::pair<int, int> A22_start = std::make_pair(A_start.first + n/2, A_start.second + n/2); // bottom right

				std::pair<int, int> B11_start = B_start; // top left
				std::pair<int, int> B12_start = std::make_pair(B_start.first, B_start.second + n/2); // top right
				std::pair<int, int> B21_start = std::make_pair(B_start.first + n/2, B_start.second); // bottom left
				std::pair<int, int> B22_start = std::make_pair(B_start.first + n/2, B_start.second + n/2); // bottom right

				std::pair<int, int> C11_start = C_start; // top left
				std::pair<int, int> C12_start = std::make_pair(C_start.first, C_start.second + n/2); // top right
				std::pair<int, int> C21_start = std::make_pair(C_start.first + n/2, C_start.second); // bottom left
				std::pair<int, int> C22_start = std::make_pair(C_start.first + n/2, C_start.second + n/2); // bottom right

				// Recursively multiply the submatrices
				_recur_matrix(A, B, C, n/2, A11_start, B11_start, C11_start);
				_recur_matrix(A, B, C, n/2, A11_start, B12_start, C12_start);

				_recur_matrix(A, B, C, n/2, A21_start, B11_start, C21_start);
				_recur_matrix(A, B, C, n/2, A21_start, B12_start, C22_start);

				_recur_matrix(A, B, C, n/2, A12_start, B21_start, C11_start);
				_recur_matrix(A, B, C, n/2, A12_start, B22_start, C12_start);

				_recur_matrix(A, B, C, n/2, A22_start, B21_start, C21_start);
				_recur_matrix(A, B, C, n/2, A22_start, B22_start, C22_start);
			
			}
		}

		// Matrix multiplication using the recursive method
		// Recursive methods need to cut the matrix into smaller submatrices of size n/2
		// This method pads the matrices to the next power of 2 and then multiplies the matrices
		// arguments: Matrix m
		// returns: result matrix of the multiplication
		// example: 
		// | 1 2 3 |
		// | 4 5 6 |
		// | 7 8 9 |
		// output:
		// | 1 2 3 0 |
		// | 4 5 6 0 |
		// | 7 8 9 0 |
		// | 0 0 0 0 |
		// The matrix is padded to the next power of 2
		Matrix pad_matrix(Matrix &m){
			int rows = m.get_rows();
			int cols = m.get_cols();
			int max_dim = std::max(rows, cols);// This work for all matrices of any size
			int n = 1;

			// Pad the matrix to the next power of 2
			while(n < max_dim){
				n *= 2;
			}

			// Create a new matrix with the padded dimensions
			Matrix padded_matrix(n, n);
			for(int i=0; i<rows; i++){
				for(int j=0; j<cols; j++){
					padded_matrix(i,j) = m(i,j);
				}
			}
			return padded_matrix;
		}
		// The requirement is to create function headers as followed
		// double* recursive_mult(int n, double* A, double* B)

		// Matrix multiplication using the recursive method
		// This is called by the multiply_matrix method
		// but outside users can also call this method
		// arguments: Matrix m
		// returns: result matrix of the multiplication
		// The matrix object A is the left matrix and the matrix object m is the right matrix
		// example: A.recursive_mult(B);

		Matrix recursive_mult(Matrix &m){
			if(this->cols != m.rows){
				std::cout << "Matrix multiplication not possible" << std::endl;
				return Matrix(0, 0);
			}

			// Pad the matrices to the next power of 2
			Matrix padded_matrix_A = pad_matrix(*this);
			Matrix padded_matrix_B = pad_matrix(m);

			// Create a result matrix with the padded dimensions
			Matrix result_matrix(padded_matrix_A.get_rows(), padded_matrix_B.get_cols());

			// Multiply the matrices using the recursive method
			_recur_matrix(padded_matrix_A, padded_matrix_B, result_matrix, padded_matrix_A.get_rows());

			// Remove padding from result matrix
			Matrix final_result_matrix(this->rows, m.cols);
			for(int i=0; i<this->rows; i++){
				for(int j=0; j<m.cols; j++){
					final_result_matrix(i,j) = result_matrix(i,j);
				}
			}
			return final_result_matrix;
		}



		private:
		// Strassen's matrix multiplication helper function
		// This function is called by the strassen_mul method
		// arguments: Matrix A, Matrix B, n
		// returns: result matrix of the multiplication

		Matrix _strassen_mul(Matrix &A, Matrix &B, int n){
			// Original way
			if(A.get_rows() == 1){  
				Matrix C(1, 1);
				C(0, 0) = A(0, 0) * B(0, 0);
				return C;

			// Optimized way for small matrices => significantly faster even than pure iterative
			// if(A.get_rows() <= 64){ 
			// 	return A.iterative_mul(B);
			}else{
				// the orignal matrix is padded to the next power of 2 before calling this function
				// so n is always a power of 2
				int half_rows = A.get_rows() / 2;

				// Divide the matrices into 4 submatrices A11, A12, A21, A22 of size n/2
				Matrix A11 = A.get_sub_matrix(0, 0, half_rows-1, half_rows-1);
				Matrix A12 = A.get_sub_matrix(0, half_rows, half_rows-1, A.get_cols()-1);
				Matrix A21 = A.get_sub_matrix(half_rows, 0, A.get_rows()-1, half_rows-1);
				Matrix A22 = A.get_sub_matrix(half_rows, half_rows, A.get_rows()-1, A.get_cols()-1);

				// Divide the matrices into 4 submatrices B11, B12, B21, B22 of size n/2
				Matrix B11 = B.get_sub_matrix(0, 0, half_rows-1, half_rows-1);
				Matrix B12 = B.get_sub_matrix(0, half_rows, half_rows-1, B.get_cols()-1);
				Matrix B21 = B.get_sub_matrix(half_rows, 0, B.get_rows()-1, half_rows-1);
				Matrix B22 = B.get_sub_matrix(half_rows, half_rows, B.get_rows()-1, B.get_cols()-1);

				// Calculate the 10 sums
				Matrix S1 = B12 - B22;
				Matrix S2 = A11 + A12;
				Matrix S3 = A21 + A22;
				Matrix S4 = B21 - B11;
				Matrix S5 = A11 + A22;
				Matrix S6 = B11 + B22;
				Matrix S7 = A12 - A22;
				Matrix S8 = B21 + B22;
				Matrix S9 = A11 - A21;
				Matrix S10 = B11 + B12;

				// Calculate the 7 products
				Matrix P1 = _strassen_mul(A11, S1, half_rows);
				Matrix P2 = _strassen_mul(S2, B22, half_rows);
				Matrix P3 = _strassen_mul(S3, B11, half_rows);
				Matrix P4 = _strassen_mul(A22, S4, half_rows);
				Matrix P5 = _strassen_mul(S5, S6, half_rows);
				Matrix P6 = _strassen_mul(S7, S8, half_rows);
				Matrix P7 = _strassen_mul(S9, S10, half_rows);

				// Calculate the 4 quadrants of the result matrix
				Matrix C11 = P5 + P4 - P2 + P6;
				Matrix C12 = P1 + P2;
				Matrix C21 = P3 + P4;
				Matrix C22 = P5 + P1 - P3 - P7;


				// Combine the 4 quadrants into the result matrix
				Matrix C(A.get_rows(), B.get_cols());
				for(int i=0; i<half_rows; i++){
					for(int j=0; j<half_rows; j++){
						C(i,j) = C11(i,j); // top left
						C(i, j+half_rows) = C12(i,j); // top right half_rows is the same as half_cols = n/2
						C(i+half_rows, j) = C21(i, j); // bottom left
						C(i+half_rows, j+half_rows) = C22(i,j); // bottom right
					}
				}

				return C;
			}
		}

		public:
		// The requirement is to create function headers as followed
		// double* strassen_mult(int n, double* A, double* B)

		// Matrix multiplication using the Strassen method
		// This is called by the multiply_matrix method
		// but outside users can also call this method
		// arguments: Matrix m
		// returns: result matrix of the multiplication
		// The matrix object A is the left matrix and the matrix object m is the right matrix
		// example: A.strassen_mul(B);
		Matrix strassen_mul(Matrix &m){
			if(this->cols != m.rows){
				std::cout << "Matrix multiplication not possible" << std::endl;
				return Matrix(0, 0);
			}
			

			Matrix padded_matrix_A = pad_matrix(*this);
			Matrix padded_matrix_B = pad_matrix(m);

			Matrix result_matrix = _strassen_mul(padded_matrix_A, padded_matrix_B, padded_matrix_A.get_rows());

			// Remove padding from result matrix
			Matrix final_result_matrix = result_matrix.get_sub_matrix(0, 0, this->rows-1, m.cols-1);
			return final_result_matrix;
		}
};

double* row_to_col_major(double* data, int rows, int cols){
	double* data_col_major = new double[rows*cols];
	for(int i=0; i<rows; i++){
		for(int j=0; j<cols; j++){
			data_col_major[j*rows + i] = data[i*cols + j];
		}
	}
	return data_col_major;
}


// Read the file and return the data
// This is dessigned for the assignment that all matrices are square
// The first line is the first matrix in column major order
// The second line is the second matrix in column major order
//
// arguments: file path
// returns: tuple of n, data1, data2
// n is the number of rows or columns in the matrix
// data1 is the first matrix in column major order
// data2 is the second matrix in column major order

std::tuple<int, double*, double*> read_file(const std::string &filename){
	
	std::ifstream file(filename, std::ios::in);
	std::string line;
	std::vector<double> data1;
	std::vector<double> data2;

	if(!file){
		std::cout << "File not found" << std::endl;
		return std::make_tuple(0, nullptr, nullptr);
	}
	//print the file content

	if(std::getline(file, line) && !line.empty()){
		std::istringstream iss(line);
		double num;
		//print line
		std::cout << line << std::endl;
		while(iss >> num){
			data1.push_back(num);
		}
	}

	// Read the second matrix
	if(std::getline(file, line) && !line.empty()){
		std::istringstream iss(line);
		//print line
		std::cout << line << std::endl;
		double num;
		while(iss >> num){
			data2.push_back(num);
		}
	}


	double* data1_arr = new double[data1.size()];
	double* data2_arr = new double[data2.size()];
	for(int i=0; i<data1.size(); i++){
		data1_arr[i] = data1[i];
	}
	for(int i=0; i<data2.size(); i++){
		data2_arr[i] = data2[i];
	}

	// return n, data1, data2
	return std::make_tuple((int)sqrt(data1.size()), data1_arr, data2_arr);
}

void test_case(){
	std::string file_1 = "Assignment 4 Files/Matrix2x2.txt";
	std::string file_2 = "Assignment 4 Files/Matrix4x4.txt";
	std::string file_3 = "Assignment 4 Files/Matrix10x10.txt";

	std::string files[] = {file_1, file_2, file_3};

	for(int i=0; i<3; i++){
		// Get the filename
		std::string filename = files[i].substr(files[i].find_last_of("/")+1);
		std::cout << "Filename: " << filename << std::endl;

		// Read the file
		std::tuple<int, double*, double*> data = read_file(files[i]);
		
		// Get the data
		int n =  std::get<0>(data);
		double* data1 = std::get<1>(data);
		double* data2 = std::get<2>(data);

		// Create the matrices
		Matrix<double> m1(n, n, data1);
		Matrix<double> m2(n, n, data2);

		// Multiply the matrices using the iterative, recursive, and Strassen methods
		double iter_time = m1.multiply_matrix(m2, "iter");
		double recur_time = m1.multiply_matrix(m2, "recur");
		double strassen_time = m1.multiply_matrix(m2, "strassen");
		
	}


}
void print_progress_bar(int current, int total, int bar_width = 50) {
    float progress = (float)current / total;
    int pos = bar_width * progress;

    std::cout << "[";
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) std::cout << "#";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}

void run_test(bool verbose=true){
	// Test the matrix multiplication on the following sizes
	int sizes[] = {2, 4, 8, 16, 32, 64, 128, 256 , 512};

	//int sizes[] = {2, 4, 6, 8, 10, 13, 15};
	int sizes_len = sizeof(sizes)/sizeof(sizes[0]);

	// Store the runtimes for each method
	std::vector<double*> runtimes;
	int longest_element = 0;
	// Loop through the sizes
	for(int i=0; i<sizes_len; i++){
		if (verbose){
			std::cout<<std::endl << "Size: " << sizes[i] << "x" << sizes[i] << std::endl << std::endl;
		}
		// Create the matrices data
		int n = sizes[i];
		double* data1 = new double[n*n];
		double* data2 = new double[n*n];

		for(int i=0; i<n*n; i++){
			data1[i] = i + 1;
			data2[i] =  n*n + i + 1;
		}

		// Create the matrices
		Matrix<double> m1(n, n, data1);
		Matrix<double> m2(n, n, data2);

		double *run_time_i = new double[3];


		
		run_time_i[0] = m1.multiply_matrix(m2, "iter", verbose);
		run_time_i[1] = m1.multiply_matrix(m2, "recur",  verbose);
		run_time_i[2] = m1.multiply_matrix(m2, "strassen", verbose);


		// Find the longest element in the matrix
		// to format the matrix printing
		for (int j=0; j<3; j++){
			std::ostringstream oss;
			oss << std::fixed << std::setprecision(6) << run_time_i[j];
			if(oss.str().length() > longest_element){
				longest_element = oss.str().length();
			}
		}
		runtimes.push_back(run_time_i);

		// Print the progress bar
		print_progress_bar(i + 1, sizes_len);

	}

	// Print the runtimes
	std::cout << std::endl;
	longest_element = longest_element * 3 / 2;
	std::cout << std::setw(longest_element) << "Size" 
	          << std::setw(longest_element) << "Iter" 
	          << std::setw(longest_element) << "Recur" 
	          << std::setw(longest_element) << "Strass" 
	          << std::endl;
	for(int i=0; i< sizes_len; i++){
		std::ostringstream oss1, oss2, oss3;
		oss1 << std::fixed << std::setprecision(8) << runtimes[i][0];
		oss2 << std::fixed << std::setprecision(8) << runtimes[i][1];
		oss3 << std::fixed << std::setprecision(8) << runtimes[i][2];
		std::cout << std::setw(longest_element) << sizes[i] 
		          << std::setw(longest_element) << oss1.str() 
		          << std::setw(longest_element) << oss2.str() 
		          << std::setw(longest_element) << oss3.str() 
		          << std::endl;
	}
	
}
int main(int argc, char const *argv[]){
	if (argc != 2) {		
		std::cout << "Usage: " << argv[0] << " <option 1/2>" << std::endl;
		std::cout << "Options: " << std::endl;
		std::cout << "1 - Test case" << std::endl;
		std::cout << "2 - Run Test on 2x2, 4x4, 8x8, 16x16, 32x32, 64x64, 128x128, 256x256, 512x512 matrices" << std::endl;
		return 1;
	}
	int option = std::stoi(argv[1]);
	if(option == 1){
		test_case();
	}else if(option == 2){
		run_test(false);
	}
}