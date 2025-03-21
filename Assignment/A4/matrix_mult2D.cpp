#include <iostream>
#include <vector>
#include <iomanip>
#include <chrono>
// use pair
#include <utility>

// Class for matrix
template <typename T>
class Matrix{
	private:
		std::vector<std::vector<T>> matrix;
		int rows;
		int cols;
	public:
		Matrix(int rows, int cols, std::vector<T> data){
			this->rows = rows;
			this->cols = cols;
			matrix.resize(rows, std::vector<T>(cols));
			for(int i=0; i<rows; i++){
				for(int j=0; j<cols; j++){
					matrix[i][j] = data[i*cols + j];
				}
			}
		}

		Matrix(int rows, int cols){
			this->rows = rows;
			this->cols = cols;

			// Initialize matrix with zeros
			matrix.resize(rows, std::vector<T>(cols));

			for(int i=0; i<rows; i++){
				for(int j=0; j<cols; j++){
					matrix[i][j] = 0;
				}
			}
		}
		
		Matrix get_sub_matrix(int start_row, int start_col, int end_row, int end_col){
			int sub_rows = end_row - start_row + 1;
			int sub_cols = end_col - start_col + 1;
			Matrix sub_matrix(sub_rows, sub_cols);

			for(int i=0; i<sub_rows; i++){
				for(int j=0; j<sub_cols; j++){
					sub_matrix[i][j] = matrix[start_row + i][start_col + j];
				}
			}
			return sub_matrix;
		}
		//Getters
		int get_rows(){
			return this->rows;
		}

		int get_cols(){
			return this->cols;
		}

		//Setters
		void set_rows(int rows){
			this->rows = rows;
		}

		void set_cols(int cols){
			this->cols = cols;
		}

		// Set operations [][] value to matrix
		std::vector<T>& operator[](int i){
			// add check for out of bounds later
			return matrix[i]; 
		}

		// Matrix addition
		Matrix operator+(Matrix &m){
			if(this->rows != m.rows || this->cols != m.cols){
				std::cout << "Matrix addition not possible" << std::endl;
				return Matrix(0, 0, std::vector<T>());
			}
			Matrix result_matrix(this->rows, this->cols);

			for(int i=0; i<this->rows; i++){
				for(int j=0; j<this->cols; j++){
					result_matrix[i][j] = this->matrix[i][j] + m.matrix[i][j];
				}
			}
			return result_matrix;
		}

		// Matrix subtraction
		Matrix operator-(Matrix &m){
			if(this->rows != m.rows || this->cols != m.cols){
				std::cout << "Matrix subtraction not possible" << std::endl;
				return Matrix(0, 0, std::vector<T>());
			}
			Matrix result_matrix(this->rows, this->cols);

			for(int i=0; i<this->rows; i++){
				for(int j=0; j<this->cols; j++){
					result_matrix[i][j] = this->matrix[i][j] - m.matrix[i][j];
				}
			}
			return result_matrix;
		}


		//Print matrix
		void print(){
			int longest_element = 0;
			for(int i=0; i<this->rows; i++){
				for(int j=0; j<this->cols; j++){
					if(std::to_string(matrix[i][j]).length() > longest_element){
						longest_element = std::to_string(matrix[i][j]).length();
					}
				}
			}

			for(int i=0; i<this->rows; i++){
				std::cout << "| ";
				for(int j=0; j<this->cols; j++){
					std::cout<< std::setw(longest_element) << matrix[i][j] << " ";
				}
				std::cout << "|";
				std::cout << std::endl;
			}
		}

		void print_matrix_multiplication(Matrix &A, Matrix &B, Matrix &result){
			// This function prints the matrices A, B and the result of their multiplication
			// Assuming result matrix is already calculated from multiplying A and B

			A.print();
			std::cout  << "*" << std::endl;
			B.print();
			std::cout  << "=" << std::endl;
			result.print();

		}


		Matrix multiply_matrix(Matrix &m, std::string method="iter"){
			if(method == "iter"){
				auto start = std::chrono::high_resolution_clock::now();
				Matrix result = this->iter_mul(m);
				auto end = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed = end - start;
				std::cout << "Time taken for matrix multiplication using iterative method: " << elapsed.count() << "s" << std::endl;

				this->print_matrix_multiplication(*this, m, result);
				return result;
			}
			else if(method == "recur"){
				auto start = std::chrono::high_resolution_clock::now();
				Matrix result = this->recur_mul(m);
				auto end = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed = end - start;
				std::cout << "Time taken for matrix multiplication using recursive method: " << elapsed.count() << "s" << std::endl;

				this->print_matrix_multiplication(*this, m, result);
				return result;
			}
			else if(method == "strassen"){
				auto start = std::chrono::high_resolution_clock::now();
				Matrix result = this->strassen_mul(m);
				auto end = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed = end - start;
				std::cout << "Time taken for matrix multiplication using strassen method: " << elapsed.count() << "s" << std::endl;

				this->print_matrix_multiplication(*this, m, result);
				return result;
			}
			else{
				std::cout << "Invalid method" << std::endl;
			}
			return Matrix(0, 0, std::vector<T>());
		}
		Matrix iter_mul(Matrix &m){
			if(this->cols != m.rows){
				std::cout << "Matrix multiplication not possible" << std::endl;
				return Matrix(0, 0, std::vector<T>());
			}
			Matrix result_matrix(this->rows, m.cols);


			for(int i=0; i<this->rows; i++){
				for(int j=0; j<m.cols; j++){
					T sum = 0;
					for(int k=0; k<this->cols; k++){
						sum += this->matrix[i][k] * m.matrix[k][j];
					}
					result_matrix.matrix[i][j] = sum;
				}
			}
			return result_matrix;
		}

		void _recur_matrix(Matrix &A, Matrix &B, Matrix &C, int n, 
			std::pair<int, int> A_start={0,0}, std::pair<int, int> B_start={0,0}, std::pair<int, int> C_start={0,0}){

			if (n == 1){
				C[C_start.first][C_start.second] += A[A_start.first][A_start.second] * B[B_start.first][B_start.second];
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
		Matrix pad_matrix(Matrix &m){
			int rows = m.get_rows();
			int cols = m.get_cols();
			int max_dim = std::max(rows, cols);
			int n = 1;
			while(n < max_dim){
				n *= 2;
			}

			Matrix padded_matrix(n, n);
			for(int i=0; i<rows; i++){
				for(int j=0; j<cols; j++){
					padded_matrix[i][j] = m[i][j];
				}
			}
			return padded_matrix;
		}
		Matrix recur_mul(Matrix &m){
			if(this->cols != m.rows){
				std::cout << "Matrix multiplication not possible" << std::endl;
				return Matrix(0, 0, std::vector<T>());
			}

			Matrix padded_matrix_A = pad_matrix(*this);
			Matrix padded_matrix_B = pad_matrix(m);

			Matrix result_matrix(padded_matrix_A.get_rows(), padded_matrix_B.get_cols());

			_recur_matrix(padded_matrix_A, padded_matrix_B, result_matrix, padded_matrix_A.get_rows());

			// Remove padding from result matrix
			Matrix final_result_matrix(this->rows, m.cols);
			for(int i=0; i<this->rows; i++){
				for(int j=0; j<m.cols; j++){
					final_result_matrix[i][j] = result_matrix[i][j];
				}
			}
			return final_result_matrix;
		}
		Matrix _strassen_mul(Matrix &A, Matrix &B, int n){
			if(A.get_rows() == 1){
				Matrix C(1, 1);
				C[0][0] = A[0][0] * B[0][0];
				return C;
			}else{
				int half_rows = A.get_rows() / 2; // Assuming A is a square matrix, A, B have same dimensions and are powers of 2
				
				// Divide the matrices into 4 submatrices
				Matrix A11 = A.get_sub_matrix(0, 0, half_rows-1, half_rows-1);
				Matrix A12 = A.get_sub_matrix(0, half_rows, half_rows-1, A.get_cols()-1);
				Matrix A21 = A.get_sub_matrix(half_rows, 0, A.get_rows()-1, half_rows-1);
				Matrix A22 = A.get_sub_matrix(half_rows, half_rows, A.get_rows()-1, A.get_cols()-1);

				// Divide the matrices into 4 submatrices
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
						C[i][j] = C11[i][j];
						C[i][j+half_rows] = C12[i][j];
						C[i+half_rows][j] = C21[i][j];
						C[i+half_rows][j+half_rows] = C22[i][j];
					}
				}

				return C;
			}
		}

		Matrix strassen_mul(Matrix &m){
			if(this->cols != m.rows){
				std::cout << "Matrix multiplication not possible" << std::endl;
				return Matrix(0, 0, std::vector<T>());
			}

			Matrix padded_matrix_A = pad_matrix(*this);
			Matrix padded_matrix_B = pad_matrix(m);

			Matrix result_matrix = _strassen_mul(padded_matrix_A, padded_matrix_B, padded_matrix_A.get_rows());

			// Remove padding from result matrix
			Matrix final_result_matrix = result_matrix.get_sub_matrix(0, 0, this->rows-1, m.cols-1);
			return final_result_matrix;
		}
};

int main(int argc, char const *argv[]){
	// std::vector<int> data1 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
	// std::vector<int> data2 = {7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22};
	std::vector<int> data1 = {
		1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
	   11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
	   21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
	   31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
	   41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
	   51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
	   61, 62, 63, 64, 65, 66, 67, 68, 69, 70,
	   71, 72, 73, 74, 75, 76, 77, 78, 79, 80,
	   81, 82, 83, 84, 85, 86, 87, 88, 89, 90,
	   91, 92, 93, 94, 95, 96, 97, 98, 99,100
   };

   std::vector<int> data2 = {
		1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
	   11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
	   21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
	   31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
	   41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
	   51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
	   61, 62, 63, 64, 65, 66, 67, 68, 69, 70,
	   71, 72, 73, 74, 75, 76, 77, 78, 79, 80,
	   81, 82, 83, 84, 85, 86, 87, 88, 89, 90,
	   91, 92, 93, 94, 95, 96, 97, 98, 99,100
   };
//    Matrix<int> m1(10, 10, data1);
//    Matrix<int> m2(10, 10, data2);

   if (argc != 2) {
		std::cout << "Usage: " << argv[0] << "<matrix_size>" << std::endl;
		return 1;
	}

   int n = std::stoi(argv[1]);
   std::vector<double> data3(n*n);
   for(int i=0; i< n*n; i++){
	   data3[i] = i+1;
   }



   Matrix<double> m1(n, n, data3);
   Matrix<double> m2(n, n, data3);

	m1.multiply_matrix(m2, "iter");
	m1.multiply_matrix(m2, "recur");
	m1.multiply_matrix(m2, "strassen");   	

	std::cout<< "Matrix size: " << n << "x" << n << std::endl;
	return 0;
}