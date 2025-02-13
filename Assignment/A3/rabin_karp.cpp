// Huy Nguyen
// CS 6420
// Assignment #3
// Dr. Rague
// Due: 02/13/2025
// Version: 1.0
// -----------------------------------------------------------------
// This program implements the Rabin-Karp 2D algorithm to search for
// a pattern matrix within a text matrix. It includes a matrix class
// to handle matrix operations and a RabinKarp2D class to perform
// the search.
//
// ANALYSIS
// -----------------------------------------------------------------
//
// This program implements the Rabin-Karp 2D algorithm based on the idea of classic Rabin-Karp, which is on 1D.
// There is the pattern matrix(MxM) and the text matrix(NxN). Given the number of characters in the alphabet d and a prime number q.
//
// 1. Pattern Hashing
// The pattern hash is done by performing hash on the pattern matrix row-wise then column-wise.
//      - step 1: row-wise hashing: for each column, it loops through each row collapsing it into matrix (1xM), each column have M elemnts so O(M). We have M columns => O(M^2)
//      - step 2: column-wise hashing: We loop the matrix 1xM collapsing it into (1x1) matrix. We have M rows => O(M)
// There is a vector declaration which takes O(M^2).
// The complexity of hashing the pattern matrix is O(M^2) + O(M) + O(M^2) = O(M^2)
//
// 2. Text Hashing	
// The text hashing is also done by performing hash on the text matrix row-wise then column-wise. But it is done on rolling hash technique.
// The idea is to hash the first window of size M of any row or column. where M < N
// Then move to the next window by removing the first element and adding the next element, which is done in arithmetic operation O(1) and accessing the element O(1).
// The step is similar:
//      - step 1: row-wise hashing: it hashes the first window of size M of each column. then it iteratively move(N-M times) to the next window by removing the first element and adding the next element.
//                the addition and removal of the element is done in O(1) and accessing the element is O(1). We have N rows => O(M + N - M) = O(N)
//                This results in matrix of size (N-M+1)xN
//      - step 2: column-wise hashing: it hashes the first window of size M of each row. then it iteratively move(N-M times) to the next window by removing the first element and adding the next element.
//                Similarly, addition and removal is constant. Because of matrix (N-M+1)xN, we have N columns but only (N-M+1) rows => O(N - M + 1) * O(N) = O(N^2), because N > M
// There are two matrix declarations which are row-wise-O(N*(N-M+1)) = O(N^2) and column-wise-O(N-M+1)*O(N - M + 1) = O(N^2)
// The complexity of hashing the text matrix is O(N^2) = O(N^2)
//
//
// 3. Searching
// Once we have hashed the pattern matrix and the text matrix, we can now search for the pattern matrix in the text matrix.
// The text matrix has size (N-M+1)x(N-M+1) and the pattern matrix is one value.
// we loop through each element of the text matrix and compare it with the hashed pattern matrix.
// the comparison is done in O(1) because we are comparing one value.
// Once we find a match, we check if the submatrix matches the pattern matrix. This is done by getting the submatrix of the text matrix and comparing it with the pattern matrix.
// Because the submatrix is of size MxM, we have to compare each element of the submatrix with the pattern matrix. This is done in O(M^2)
// However, the hashing is done in the manner that reduces the number of collisions, which reduces the number of false positives.
// So we can assumes that the sub matrix comparison is done, in most cases, only once O(MxM)
// So it is O((N-M+1)^2) - 1 + O(M^2) = O(N^2)
// The complexity of searching is O(N^2)


// -----------------------------------------------------------------

// Compiler directives
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <utility>
#include <random>
#include <chrono>
#include <fstream>
#include <climits>
#include <iomanip> // For std::setw


// -----------------------------------------------------------------
// This class represents a matrix and provides various operations
// to manipulate and retrieve data from the matrix.
//
// Version 1.0
// ----------------------------------------------------------

// Template-based matrix class
template <typename T>
class matrix {
private:
	// 2D vector to store matrix data
	std::vector<std::vector<T>> data;
	// Dimensions of the matrix
	int n1, n2;

public:
	// Default constructor
	matrix() : n1(0), n2(0) {}

	// Constructor to initialize matrix with given data
	matrix(std::vector<std::vector<T>> data) {
		this->data = data;
		this->n1 = data.size();
		this->n2 = data[0].size();
	}

	// Constructor to initialize matrix with given dimensions and fill with zeros
	matrix(int n1, int n2) {
		this->n1 = n1;
		this->n2 = n2;
		this->data = std::vector<std::vector<T>>(n1, std::vector<T>(n2, 0));
	}

	// Set matrix data and update dimensions
	void set_data(std::vector<std::vector<T>> data) {
		this->data = data;
		this->n1 = data.size();
		this->n2 = data[0].size();
	}

	// Get matrix data
	std::vector<std::vector<T>> get_data() {
		return this->data;
	}

	// Get matrix dimensions
	std::pair<int, int> get_size() {
		return std::pair<int, int>(this->n1, this->n2);
	}

	// Get element at position (i, j)
	T get_element(int i, int j) {
		return this->data[i][j];
	}

	// Set element at position (i, j)
	void set_element(int i, int j, T value) {
		this->data[i][j] = value;
	}

	// Set entire row at index i
	void set_row(int i, std::vector<T> row) {
		this->data[i] = row;
	}

	// Set entire column at index j
	void set_col(int j, std::vector<T> col) {
		for (int i = 0; i < this->n1; i++) {
			this->data[i][j] = col[i];
		}
	}

	// Get entire row at index i
	std::vector<T> get_row(int i) {
		return this->data[i];
	}

	// Get entire column at index j
	std::vector<T> get_col(int j) {
		std::vector<T> col;
		for (int i = 0; i < this->n1; i++) {
			col.push_back(this->data[i][j]);
		}
		return col;
	}

	// Get submatrix from (i1, j1) to (i2-1, j2-1)
	matrix<T> get_submatrix(int i1, int j1, int i2, int j2) {
		// Check for valid submatrix indices
		if (i1 < 0 || i2 > this->n1 || j1 < 0 || j2 > this->n2) {
			std::cout << "Invalid submatrix" << std::endl;
			return matrix<T>();
		}

		std::vector<std::vector<T>> submatrix;
		for (int i = i1; i < i2; i++) {
			std::vector<T> row;
			for (int j = j1; j < j2; j++) {
				row.push_back(this->data[i][j]);
			}
			submatrix.push_back(row);
		}

		return matrix<T>(submatrix);
	}

	// Get a random submatrix with given dimensions
	matrix<T> get_submatrix_random(int dim1, int dim2) {
		// Check for valid submatrix indices
		if (dim1 > this->n1 || dim2 > this->n2) {
			std::cout << "Invalid submatrix" << std::endl;
			return matrix<T>();
		}

		// Get time seed at now
		srand(time(0));

		// Get random indices that x + dim1 < n1 and y + dim2 < n2 and x, y > 0
		int x = (this->n1 - dim1) > 0 ? rand() % (this->n1 - dim1) : 0;
		int y = (this->n2 - dim2) > 0 ? rand() % (this->n2 - dim2) : 0;

		return get_submatrix(x, y, x + dim1, y + dim2);
	}

	// Get submatrix as a vector of vectors
	std::vector<std::vector<T>> get_submatrix_vector(int i1, int j1, int i2, int j2) {
		// Check for valid submatrix indices
		if (i1 < 0 || i2 > this->n1 || j1 < 0 || j2 > this->n2) {
			std::cout << "Invalid submatrix" << std::endl;
			return std::vector<std::vector<T>>();
		}

		std::vector<std::vector<T>> submatrix;
		for (int i = i1; i < i2; i++) {
			std::vector<T> row;
			for (int j = j1; j < j2; j++) {
				row.push_back(this->data[i][j]);
			}
			submatrix.push_back(row);
		}

		return submatrix;
	}

	// Get a random submatrix as a vector of vectors
	std::vector<std::vector<T>> get_random_submatrix(int dim1, int dim2) {
		// Check for valid submatrix indices
		if (dim1 > this->n1 || dim2 > this->n2) {
			std::cout << "Invalid submatrix" << std::endl;
			return std::vector<std::vector<T>>();
		}

		std::vector<std::vector<T>> submatrix;

		// Get time seed at now
		srand(time(0));

		// Get random indices that x + dim1 < n1 and y + dim2 < n2 and x, y > 0
		int x = (this->n1 - dim1) > 0 ? rand() % (this->n1 - dim1) : 0;
		int y = (this->n2 - dim2) > 0 ? rand() % (this->n2 - dim2) : 0;

		for (int i = x; i < x + dim1; i++) {
			std::vector<T> row;
			for (int j = y; j < y + dim2; j++) {
				row.push_back(this->data[i][j]);
			}
			submatrix.push_back(row);
		}

		return submatrix;
	}

	// Print matrix to console
	void print(std::string sep = "") {
		if (!sep.empty()) {
			int max_width = 0;

			// Calculate the maximum width of the elements
			for (int i = 0; i < this->n1; i++) {
				for (int j = 0; j < this->n2; j++) {
					int width = std::to_string(this->data[i][j]).length();
					if (width > max_width) {
						max_width = width;
					}
				}
			}

			// Print the matrix with proper alignment
			for (int i = 0; i < this->n1; i++) {
				for (int j = 0; j < this->n2; j++) {
					std::cout << std::left << std::setw(max_width) << this->data[i][j];
					if (j < this->n2 - 1) {
						std::cout << sep;
					}
				}
				std::cout << std::endl;
			}
		} else {
			// Print the matrix without alignment
			for (int i = 0; i < this->n1; i++) {
				for (int j = 0; j < this->n2; j++) {
					std::cout << this->data[i][j] << sep;
				}
				std::cout << std::endl;
			}
		}
	}
};


// -----------------------------------------------------------------
// This class represents the Rabin-Karp 2D algorithm and provides
// various operations to perform the search for a pattern matrix
// within a text matrix.
//
// Version 1.0
// ----------------------------------------------------------
class RabinKarp2D {
private:
    // Declare matrix
    matrix<long> T; // Text matrix
    matrix<long> P; // Pattern matrix

    matrix<long> hashed_T; // Hashed text matrix
    long hashed_P; // Hashed pattern matrix

    int d; // Number of characters in the alphabet
    int q; // Prime number
	int h; // d^(m -1) 
public:
    // Function to return a random large prime number
    int random_large_prime() {
        int min = 0, max = 49;
        long primes[] = {
            2094665479L, 1783990163L, 2094521287L, 2134397081L, 2126326253L, 
            1957216747L, 1436547389L, 1428780767L, 2075625529L, 1593123733L, 
            2132587157L, 1965562429L, 1164701777L, 1568991883L, 2130061793L, 
            1075370311L, 1711832929L, 2054631589L, 1587361861L, 1435348609L, 
            1332084959L, 1465215911L, 2088173753L, 1933073123L, 1319415599L, 
            1211741129L, 1487473783L, 1656920599L, 1817614213L, 1838911937L, 
            1697951429L, 1673793083L, 1971101663L, 1570547117L, 1869368041L, 
            1855484017L, 2057695543L, 1806695647L, 2082498797L, 2090345119L, 
            1349212999L, 1456810283L, 1271362889L, 1959057733L, 1073964823L, 
            1315871351L, 1308843649L, 1543027127L, 1230659387L, 1828780297L 
        };

        static std::default_random_engine engine;
        engine.seed(std::chrono::system_clock::now().time_since_epoch().count());
        static std::uniform_int_distribution<int> dist(min, max);

        return (int)primes[dist(engine)];
    }

    // Default constructor
    RabinKarp2D() {
        this->d = 256;
        this->q = random_large_prime();
    }

    // Constructor to initialize RabinKarp2D with given data
    RabinKarp2D(matrix<long> T, matrix<long> P, int d, int q=0) {
        this->T = T;
        this->P = P;
        this->d = d;
        this->q = (q == 0) ? random_large_prime() : q;
        this->h = (int)pow(d, P.get_size().first * P.get_size().second - 1) % q;

        // Perform hashing
        perform_hashing();
    }

    // Constructor to initialize RabinKarp2D with given dimensions
    RabinKarp2D(int d, int q): d(d), q(random_large_prime()) {
        this->d = d;
        this->q = q;
    }

private:
    // Setters
    void setT(std::vector<std::vector<long>> T) {
        this->T = T;
    }

    void setP(std::vector<std::vector<long>> P) {
        this->P = P;
    }

    void setD(int d) {
        this->d = d;
    }

    void setQ(int q) {
        this->q = q;
    }

    // Perform hashing on the given matrix
    void perform_hashing() {
        std::pair<int, int> P_size = this->P.get_size();
        int m1 = P_size.first;
        int m2 = P_size.second;

        // Hash the pattern matrix
        matrix<long> H;
        long hash_P = hash_pattern(P);
        this->hashed_P = hash_P;

        // Hash the text matrix
        // The hash_pattern is designed to run row-wise then column-wise
        H = hash(T, 0);
        H = hash(H, 1);

        // Store the hashed text matrix
        this->hashed_T = H;
    }

public:
	// Set data for the RabinKarp2D algorithm
	// input being matrix of long
    void set_data(matrix<long> T, matrix<long> P, int d, int q=0) {
        this->T = T;
        this->P = P;

        setD(d);

        // Set a random large prime number if q is not provided
        if (q == 0) {
            setQ(random_large_prime());
        } else {
            setQ(q);
        }

        // Perform hashing
        perform_hashing();
    }



	// Calculate the hash value for the first window
	// input being vector of vectors
    void set_data(std::vector<std::vector<long>> T, std::vector<std::vector<long>> P, int d, int q=0) {
        setT(T);
        setP(P);
        setD(d);

        // Set a random large prime number if q is not provided
        if (q == 0) {
            setQ(random_large_prime());
        } else {
            setQ(q);
        }

        // Perform hashing
        perform_hashing();
    }


	// Print the text and pattern matrices
    void print_matrix() {
        std::cout << "Text Matrix: " << std::endl;
        T.print();
        std::cout << std::endl;
        std::cout << "Pattern Matrix: " << std::endl;
        P.print();
        std::cout << std::endl;
    }

	// Print the hashed text and pattern matrices
    void print_hashed_matrix() {
        std::cout << "Hashed Text Matrix: " << std::endl;
        hashed_T.print(" ");
        std::cout << std::endl;
        std::cout << "Hashed Pattern Matrix: " << hashed_P << std::endl;
        std::cout << std::endl;
    }

    // Getters
    matrix<long> getT() {
        return this->T;
    }

    matrix<long> getP() {
        return this->P;
    }

    int getD() {
        return this->d;
    }

    int getQ() {
        return this->q;
    }



	////////////////////////////////////////////////////////////////
	// Hashing functions
	////////////////////////////////////////////////////////////////

    // Calculate the hash value for the first window
    long first_window_hash(matrix<long> M, int index, int axis=0) {
        // Assuming the matrix T is set
        // d, q, h are also set
        long hash = 0;

        // Get the window size either row-wise or column-wise
        int window_size = (axis == 0) ? P.get_size().first : P.get_size().second;

        // Calculate the hash value for the first window
        for (int i = 0; i < window_size; i++) {
            // Get the element from the matrix T either row-wise or column-wise
            // If axis 0, keep col move along the rows
            int element = (axis == 0) ? M.get_element(i, index) : M.get_element(index, i);
            hash = (d * hash + element) % q;
        }

        return hash;
    }

    // Hash the pattern matrix
	// it calculates the hash value of a given matrix
	// row wise then column wise 
    long hash_pattern(matrix<long> P) {
        std::vector<long> hash_values; // 

        for (int i = 0; i < P.get_size().second; i++) {
            long hash = first_window_hash(P, i, 0);
            hash_values.push_back(hash);
        }

        long overall_hash = 0;

        for (int i = 0; i < hash_values.size(); i++) {
            overall_hash = (d * overall_hash + hash_values[i]) % q;
        }

        return overall_hash;
    }


    // Assign value to the hashed matrix
    void assign_value(matrix<long> &H, int hash, int idx1, int idx2, int axis=0) {
        // Axis = 0 for row-wise, axis = 1 for column-wise
        if (axis == 0) {
            H.set_element(idx1, idx2, hash); // Row-wise
        } else {
            H.set_element(idx2, idx1, hash);
        }
    }

    // Perform Rabin-Karp rolling hash on the given matrix on the specified axis
    matrix<long> hash(matrix<long> M, int axis=0) {
        // Assuming the matrix T is set
        // d, q, h are also set
        // Axis = 0 for collapse row-wise, mean take a column and move along the rows
        // Axis = 1 for collapse column-wise mean take a row and move along the columns
        int new_dim1, new_dim2;
        new_dim1 = (axis == 0) ? M.get_size().first - P.get_size().first + 1 : M.get_size().first;
        new_dim2 = (axis == 0) ? M.get_size().second : M.get_size().second - P.get_size().second + 1;

        matrix<long> H(new_dim1, new_dim2); // Hashed matrix O(N^2)

        // If axis is 0, we are hashing row-wise meaning we take a column and move along the rows
        int first_axis = (axis == 0) ? M.get_size().second : M.get_size().first;
        int second_axis = (axis == 0) ? M.get_size().first : M.get_size().second;

        // But to get the window size, we need to know the size of the pattern matrix,
        // so axis 0 means we are hashing row-wise, we need to know the number of rows in the pattern matrix
        int first_window_size = (axis == 0) ? P.get_size().second : P.get_size().first;
        int second_window_size = (axis == 0) ? P.get_size().first : P.get_size().second;

        // Move along the second axis because if axis is 0,
        // we are hashing row-wise mean we first move along the columns
        // and take each value on the row.
        // If axis is 1, we are hashing column-wise, mean we first move along the rows
        // and take each value on the column
        for (int i = 0; i < first_axis; i++) {  // Loop through the first axis up to N O(N)
            long hash = first_window_hash(M, i, axis); // Hash the first window of size M O(M)
            long h = (int)pow(d, second_window_size - 1) % q;
            h = 1;
            for (int k = 1; k <= second_window_size - 1; k++) {
                h = (d * h) % q;
            }

            assign_value(H, hash, 0, i, axis);

            for (int j = 1; j < second_axis - second_window_size + 1; j++) { // Loop through the second axis up to n - m  O(N - M)
                // Get the first and last element of the window
                // Axis = 0 mean we are hashing row-wise
                // We freeze the column and move along the rows
                int first_element = (axis == 0) ? M.get_element(j - 1, i) : M.get_element(i, j - 1);
                // We freeze the first axis
                // And move along the second axis at j + m - 1 for the last element 1 2 3 4 5. j = 1 m = 3 => last element = 4 at index 3
                int last_element = (axis == 0) ? M.get_element(j + second_window_size - 1, i) : M.get_element(i, j + second_window_size - 1);

                hash = (hash + q - h * first_element % q) % q;
                hash = (hash * d + last_element) % q;

                assign_value(H, hash, j, i, axis); // Assign the hash value to the hashed matrix O(1)
            }
        }
        return H;
    }

    // Search for the pattern matrix within the text matrix
	// This happens after the hashing is done either by the constructor or by the set_data method
    std::pair<int, int> search() {
        for (int row = 0; row < hashed_T.get_size().first; row++) {
            for (int col = 0; col < hashed_T.get_size().second; col++) {
                if (hashed_T.get_element(row, col) == hashed_P) {
                    // Check if the submatrix matches the pattern matrix
                    matrix submatrix = T.get_submatrix(row, col, row + P.get_size().first, col + P.get_size().second);

                    if (submatrix.get_data() == P.get_data()) {
                        std::cout << "Pattern found at: " << row << ", " << col << std::endl << std::endl;
                        return std::pair<int, int>(row, col);
                    }
                }
            }
        }

        std::cout << "Pattern not found" << std::endl << std::endl;

        return std::pair<int, int>(-1, -1);
    }
};


// Function to load a matrix from a file
// The matrix should be in the format:
// 12345
// 67890
// 54321
// 09876
std::vector<std::vector<long>> load_matrix_from_file(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<std::vector<long>> matrix;
    std::string line;

    if (file.is_open()) {
        while (std::getline(file, line)) {
            std::vector<long> row;
            for (char c : line) {
                if (isdigit(c)) {
                    row.push_back(c - '0'); // Convert char to int
                }
            }
            if (!row.empty()) {
                matrix.push_back(row);
            }
        }
        file.close();
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }

    return matrix;
}




using vector2D= std::vector<std::vector<long>>;

class examples{
public:

	//constructor
	examples(){};

	//destructor
	~examples(){};

	static std::vector<std::pair<vector2D, vector2D> > get_all_examples(){
		std::vector<std::pair<vector2D, vector2D> > all_examples;
		all_examples.push_back(example_1());
		all_examples.push_back(example_2());
		all_examples.push_back(example_3());
		all_examples.push_back(example_4());
		all_examples.push_back(example_5_no_found());

		return all_examples;
	}

	//methods
	static std::pair<vector2D, vector2D> example_1(){
		vector2D T_data = {
			{4, 5, 4, 5, 1},
			{4, 0, 3, 2, 2},
			{3, 4, 5, 5, 4},
			{1, 1, 3, 0, 2},
			{2, 1, 3, 3, 2}
		};

		vector2D P_data = {
			{3, 2},
			{5, 5}
		};

		return std::pair<vector2D, vector2D>(T_data, P_data);
	}

	static std::pair<vector2D, vector2D> example_2(){
		vector2D T_data = {
			{2, 3, 9, 1, 3},
			{5, 1, 2, 3, 4},
			{2, 3, 1, 5, 3},
			{0, 5, 4, 2, 3},
			{1, 0, 6, 1, 5}
		};

		vector2D P_data = {
			{2, 3},
			{1, 5}
		};

		return std::pair<vector2D, vector2D>(T_data, P_data);
	}

	static std::pair<vector2D, vector2D> example_3() {
		vector2D T_data = {
			{4, 5, 4, 5, 1},
			{4, 0, 3, 2, 2},
			{3, 4, 5, 5, 4},
			{1, 1, 3, 0, 2},
			{2, 1, 3, 3, 2}
		};

		vector2D P_data = {
			{5, 5},
			{3, 0}
		};

		return std::pair<vector2D, vector2D>(T_data, P_data);
	}

	static std::pair<vector2D, vector2D> example_4() {
		vector2D T_data = {
			{2, 0, 5, 6, 2, 5, 2, 3, 0, 1},
			{2, 5, 5, 2, 3, 5, 3, 2, 6, 0},
			{5, 3, 5, 6, 3, 1, 5, 6, 2, 4},
			{1, 0, 2, 1, 6, 3, 0, 0, 3, 6},
			{4, 0, 2, 5, 2, 1, 0, 3, 3, 5},
			{2, 4, 0, 2, 1, 3, 5, 3, 1, 3},
			{3, 6, 1, 5, 1, 5, 6, 5, 0, 0},
			{1, 4, 0, 6, 1, 2, 4, 1, 3, 4},
			{2, 3, 1, 4, 0, 1, 1, 4, 1, 0},
			{4, 5, 0, 1, 2, 1, 2, 1, 3, 0}
		};

		vector2D P_data = {
			{3, 5, 3, 1},
			{5, 6, 5, 0},
			{2, 4, 1, 3},
			{1, 1, 4, 1}
		};

		return std::pair<vector2D, vector2D>(T_data, P_data);
	}

	static std::pair<vector2D, vector2D> example_5_no_found() {
		vector2D T_data = {
			{5, 6, 0, 3, 0},
			{3, 0, 6, 6, 0},
			{0, 0, 6, 5, 5},
			{2, 4, 3, 6, 2},
			{3, 6, 4, 6, 2}
		};

		vector2D P_data = {
			{2, 0},
			{0, 2}
		};

		return std::pair<vector2D, vector2D>(T_data, P_data);
	}
};

// -----------------------------------------------------------------
// This function runs the Rabin-Karp 2D algorithm on preset matrices
// and prints the results.
//
// Version 1.0
// ----------------------------------------------------------
void part_1_run_preset_matrices() {
    // Get the first example
    std::pair<vector2D, vector2D> example = examples::example_1();

    // Create RabinKarp2D object with given data
    RabinKarp2D rk;

    // Get all examples
    std::vector<std::pair<vector2D, vector2D>> all_examples = examples::get_all_examples();

    // Iterate through all examples and run the Rabin-Karp 2D algorithm
    for (int i = 0; i < all_examples.size(); i++) {
        std::pair<vector2D, vector2D> example = all_examples[i];
        std::cout << "=========== test " << i + 1 << " =============" << std::endl;
        rk.set_data(example.first, example.second, 256);
        rk.print_matrix();
        rk.print_hashed_matrix();
        rk.search();
        std::cout << std::endl;
    }

    std::cout << std::endl;
}

// -----------------------------------------------------------------
// This function runs the Rabin-Karp 2D algorithm on random submatrices
// and prints the results.
//
// Version 1.0
// ----------------------------------------------------------
void run_each_random_step(RabinKarp2D &rk, matrix<long> &M, int d, int t1, int t2, int p1, int p2, int print = 0) {
    // Check for valid submatrix dimensions
    if (p1 > t1 || p2 > t2) {
        std::cout << "Invalid submatrix" << std::endl;
        return;
    }

    // Get a random submatrix from the text matrix
    matrix<long> T = M.get_submatrix_random(t1, t2);

    // Get a random submatrix from within the text matrix
    matrix<long> P = T.get_submatrix_random(p1, p2);

    // Set data for the Rabin-Karp 2D algorithm
    rk.set_data(T, P, d);

    // Print the matrices if the print flag is set
    if (print) {
        rk.print_matrix();
    } else {
        std::cout << "Pattern Matrix: " << std::endl;
        rk.getP().print();
    }

    // Perform the search
    rk.search();
}

// -----------------------------------------------------------------
// This function runs the Rabin-Karp 2D algorithm on a 1000x1000 matrix
// with various submatrix sizes and prints the results.
//
// Version 1.0
// ----------------------------------------------------------
void part_2_run_on_1000_1000() {
    // Load the 1000x1000 matrix from the file
    std::vector<std::vector<long>> M_1000_1000 = load_matrix_from_file("octal.txt");
	// 	a random 2 x 2 pattern within a random 5 x 5 text.
	// 	a random 8 x 8 pattern within a random 25 x 25 text.
	// 	a random 16 x 16 pattern within a random 100 x 100 text.
	// 	a random 64 x 64 pattern within the complete 1000 x 1000 octal.txt file.
    // Create RabinKarp2D object with given data
    RabinKarp2D rk;

    // Create a matrix object from the loaded data
    matrix<long> M(M_1000_1000);

    int d = 10; // Number of characters in the alphabet

    // Run tests with various submatrix sizes
    std::cout << "=========== test 1: 2 x 2 pattern within a random 5 x 5 text =============" << std::endl;
    run_each_random_step(rk, M, d, 5, 5, 2, 2, 1);

    std::cout << "=========== test 2: 8 x 8 pattern within a random 25 x 25 text =============" << std::endl;
    run_each_random_step(rk, M, d, 25, 25, 8, 8, 1);

    std::cout << "=========== test 3: 16 x 16 pattern within a random 100 x 100 text =============" << std::endl;
    run_each_random_step(rk, M, d, 100, 100, 16, 16, 1);

    std::cout << "=========== test 4: 64 x 64 pattern within the complete 1000 x 1000 octal.txt file =============" << std::endl;
    run_each_random_step(rk, M, d, 1000, 1000, 64, 64, 0);
}
int main() {
	int d = 10; // Number of characters in the alphabet

	// print header
	std::cout << "======================================" << std::endl;
	std::cout << "Rabin-Karp 2D Part 1: Preset Matrices" << std::endl;
	std::cout << "======================================" << std::endl;
	std::cout << std::endl;
	// run test part 1
	part_1_run_preset_matrices();

	

	std::cout << "======================================" << std::endl;
	std::cout << "Rabin-Karp 2D Part 2: Random Matrices" << std::endl;
	std::cout << "======================================" << std::endl;
	std::cout << std::endl;
	// run test part 2
	part_2_run_on_1000_1000();



}
