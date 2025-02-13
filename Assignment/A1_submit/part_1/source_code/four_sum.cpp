// Huy Nguyen
// CS 6420
// Assignment #1
// Dr. Rague
// Due: 01/16/2025
// Version: 1.0
// -----------------------------------------------------------------
// SUMMARY of PART 1
// This program implements two approaches for solving the "4-Sum" problem,
// which finds all unique quadruplets of integers in an array that sum to zero.
// The program compares a brute force method with an optimized approach
// and measures their performance for different input sizes.
// -----------------------------------------------------------------
// To compile the program, use the following command:
// make
// 
// The program is designed to auto load all provided input files and run the test on them.
// To run the program, use the following command:
// ./four_sum

//-----------------------------------------------------------------
// ANALYSIS of PART 1
//
// # Brute Force Approach
// The brute force approach iterates over all possible quadruplets in the array. There are four nested loops to generate all possible combinations.
// The time complexity of this approach is O(n^4) because there are four nested loops that iterate over the array of size n.
// There are n number of unique elements in the input array, there should be nC4 = n! / (4! * (n-4)!)  = n*(n-1)*(n-2)*(n-3)/4! ~ n^4/24 possible combinations.
// Therefore, the overall time complexity of the brute force approach is O(n^4).
// Take any brute-force test run, it can be seen that the running time at 2000(x_1 * e+7) is approximately 10^4 greater than that at 200(x_2 * e+11)
// or at any running time at 1600 would be approximate 16 times greater than that at 800, 1600/800 = 2, 2^4 = 16.
//
// # Optimized Approach
// The optimized approach precomputes pairs of sums and stores them in a hash map. It then iterates through all sums and finds complementary pairs.
// There are two main loops in the optimized approach, one to generate pairs and another to find complementary pairs.
// For the first loop, because it is pairs of elements, the time complexity of this approach is O(n^2) to generate all possible pairs, which is nC2 = n*(n-1)/2 possible combinations.
// For the second loop, it iterates through all sums and finds complementary pairs. 
// It depends on the nature of the data to determine the time complexity of the second loop.
// Two of the extreme case both happens when there is a minimum number of sum-value keys(2n-3).  [0,1,2,3,4,5]
// The best case is when there is no complementary sum-value keys, and the worst case is when half of the keys are complementary to the other half [-2,-1,0,1,2]
// 
// The best case in term of complexity is when there is no complementaty sum-value keys and there is a minimum number of unique sum-value key.
// In this case, there are 2n-3 keys and they are all positive(or negative), making it O(n). 
// The whole algorithm is O(n^2) + O(n) = (n^2)
// 
// The worst case happens when there is also a minimum number of sum-value key, but in this case half of the key range should be complementary to the other half.
// This makes the second loop go through all keys (2n-3), and all pairs (nC2).
// It means that the second loop's complexity and the whole algorithm's is O(n^3)
//
// # Performance Comparison
// For the 6 provided input files, the optimized approach outperforms the brute force approach in terms of execution time.
// It is clear that the brute-force grows at a much faster rate than the optimized approach.
// Observation between 50-500, 100-1000 shows that the brute-force running time increases in a quadratic manner,
// while the optimized approach approximately O(n^2) at 50-500 and O(n^3) at 100-1000.
//
// # Detailed Analysis
// can be found in the analysis.pdf file.












#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <set>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <map>
#include <array>
#include <list>
#include "alg_stopwatch.h" // A custom stopwatch class for measuring performance
#include <thread> // Include threading library
#include <future> // Include future library for getting return values
#include <random>
#include <limits>


using namespace std;

// -----------------------------------------------------------------
// Function: four_sum_brute
// This function uses a brute force approach to solve the "4-Sum" problem.
// It iterates over all possible quadruplets in the array, checking if
// their sum is zero. It stores unique quadruplets in a set.
//
// Version: 1.0
// -----------------------------------------------------------------
int four_sum_brute(const vector<int>& arr) {
    int total_pairs = 0; // Counter for the total number of pairs found
    set<vector<int>> pairs; // A set to store unique quadruplets

    // Nested loops to iterate through all possible quadruplets
    for (int i = 0; i < arr.size(); i++) {
        for (int j = i + 1; j < arr.size(); j++) {
            for (int k = j + 1; k < arr.size(); k++) {
                for (int l = k + 1; l < arr.size(); l++) {
                    if (arr[i] + arr[j] + arr[k] + arr[l] == 0) {
                        vector<int> temp = {arr[i], arr[j], arr[k], arr[l]};
                        sort(temp.begin(), temp.end()); // Sort the quadruplet to ensure uniqueness
                        pairs.insert(temp); // Insert into the set
                        total_pairs++; // Increment the counter
                    }
                }
            }
        }
    }

    return total_pairs; // Return the count of unique quadruplets
}

// -----------------------------------------------------------------
// Function: sort_four
// This helper function sorts four integers in non-decreasing order.
// 
// Version: 1.0
// -----------------------------------------------------------------
void sort_four(int& a, int& b, int& c, int& d) {
    // Pairwise comparisons to sort the four integers
    if (a > b) swap(a, b);
    if (c > d) swap(c, d);
    if (a > c) swap(a, c);
    if (b > d) swap(b, d);
    if (b > c) swap(b, c);
}

// -----------------------------------------------------------------
// Struct: ArrayHash
// This structure provides a custom hash function for arrays of size 4.
// It is used to store unique quadruplets in an unordered_map.
// 
// Version: 1.0
// -----------------------------------------------------------------
struct ArrayHash {
    template <typename T, std::size_t N>
    std::size_t operator()(const std::array<T, N>& arr) const {
        std::size_t hash = 0; // Initialize hash value
        for (const auto& elem : arr) {
            // Combine the hash of each element using a hash combiner formula
            hash ^= std::hash<T>{}(elem) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

// -----------------------------------------------------------------
// Function: four_sum
// This function uses an optimized approach to solve the "4-Sum" problem.
// It precomputes pairs of sums and stores them in a hash map, allowing
// for faster lookups of complementary pairs. Unique quadruplets are stored
// using a custom hash function.
// 
// Version: 1.0
// -----------------------------------------------------------------
int four_sum(const vector<int>& arr) {
    unordered_map<int, vector<pair<int, int>>> sum_pairs; // Map to store pairs by their sum
    int total_pairs = 0; // Counter for the total number of pairs processed

    // Step 1: Precompute all pairs and group them by their sum
    for (int i = 0; i < arr.size(); i++) {
        for (int j = i + 1; j < arr.size(); j++) {
            sum_pairs[arr[i] + arr[j]].push_back({arr[i], arr[j]}); // Add pair to the map
            total_pairs++; // Increment the pair counter
        }
    }

    unordered_map<array<int, 4>, int, ArrayHash> quad_map; // Map to store unique quadruplets
    quad_map.reserve(total_pairs); // Reserve space to improve efficiency

    // Step 2: Iterate through all sums and find complementary pairs
    for (const auto& [sum, pairs] : sum_pairs) {
        if (sum_pairs.find(-sum) != sum_pairs.end()) { // Check if the complementary sum exists
            for (const auto& [a, b] : pairs) {
                for (const auto& [c, d] : sum_pairs[-sum]) {
                    // Ensure that no indices overlap
                    if (a != c && a != d && b != c && b != d) {
                        int a_sorted = a, b_sorted = b, c_sorted = c, d_sorted = d;
                        sort_four(a_sorted, b_sorted, c_sorted, d_sorted); // Sort to ensure uniqueness
                        array<int, 4> temp = {a_sorted, b_sorted, c_sorted, d_sorted};
                        quad_map[temp] = 1; // Store the unique quadruplet
                    }
                }
            }
        }
    }

    return quad_map.size(); // Return the count of unique quadruplets
}

// -----------------------------------------------------------------
// Function: read_file
// This function reads integers from a file and returns them as a vector.
// The file format determines whether the first line contains all values
// or if each subsequent line contains a single value.
// 
// Version: 1.0
// -----------------------------------------------------------------
vector<int> read_file(const string& filename) {
    ifstream file(filename); // Open the file
    vector<int> arr;

    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return arr;
    }

    string line;
    if (getline(file, line)) { // Read the first line
        istringstream iss(line);
        int num;
        vector<int> first_line_values;

        // Parse integers in the first line
        while (iss >> num) {
            first_line_values.push_back(num);
        }

        if (first_line_values.size() > 1) {
            arr.insert(arr.end(), first_line_values.begin(), first_line_values.end()); // Add all values
        } else if (!first_line_values.empty()) {
            arr.push_back(first_line_values[0]); // Add the single value

            // Read subsequent lines for additional values
            while (getline(file, line)) {
                istringstream iss(line);
                if (iss >> num) {
                    arr.push_back(num);
                }
            }
        }
    }

    file.close();
    return arr;
}

// -----------------------------------------------------------------
// Function: create_list_of_numbers
// This function generates a vector of random integers in the range [-n, n].
// 
// Version: 1.0
// -----------------------------------------------------------------
vector<int> create_list_of_numbers(int n) {
    vector<int> arr; // Vector to store random numbers
    srand(static_cast<unsigned>(time(0))); // Seed the random number generator

    // Generate random numbers ensuring no duplicates
    for (int i = 0; i < n; ++i) {
        int random_value = rand() % (4 * n + 1) - 2*n; // Generate number in range [-2n, 2n]
        if (find(arr.begin(), arr.end(), random_value) != arr.end()) { // Check for duplicates
            --i;
            continue;
        }
        arr.push_back(random_value);
    }

    return arr;
}
vector<int> create_list_of_number_arithmetic(int n){
	vector<int> arr; // Vector to store random numbers
    srand(static_cast<unsigned>(time(0))); // Seed the random number generator

	for (int i = 0; i < n; ++i) {
		int random_value = i; // Generate number in range [0, n-1]
		arr.push_back(random_value);
	}

	return arr;


}

std::vector<int> create_list_of_numbers_truly_random(int n) {
    std::vector<int> arr; // Vector to store unique random numbers

    // Use <random> for better random number generation
    std::random_device rd;               // Seed generator
    std::mt19937 gen(rd());              // Mersenne Twister RNG
    std::uniform_int_distribution<int> dist(std::numeric_limits<int>::min(), std::numeric_limits<int>::max());

    // Generate random numbers ensuring no duplicates
    while (arr.size() < static_cast<size_t>(n)) {
        int random_value = dist(gen);    // Generate a random number in [INT_MIN, INT_MAX]
        if (std::find(arr.begin(), arr.end(), random_value) == arr.end()) {
            arr.push_back(random_value); // Add unique random value
        }
    }

    return arr;
}

// -----------------------------------------------------------------
// Function: run_tests
// This function measures the execution time of both the brute force
// and optimized "4-Sum" solutions for different input sizes.
// 
// Version: 1.0
// -----------------------------------------------------------------
void run_tests(int start, int end, int step) {
    cout << left << setw(10) << "InputSize" << setw(20) << "OptimizedTime(ns)" << setw(20) << "BruteForceTime(ns)" << setw(20) << "OptimizedCount" << setw(20) << "BruteForceCount" << endl;	
    cout << string(90, '-') << endl;
    StopWatch watch; // Timer object for measuring execution time

    for (int size = start; size <= end; size += step) {
        vector<int> arr = create_list_of_numbers(size); // Generate input data

        // Measure execution time for the optimized solution
        watch.reset();
        int optimized_count = four_sum(arr);
        auto optimized_time = watch.elapsed_time();

        // Measure execution time for the brute force solution
        watch.reset();
        int brute_force_count = four_sum_brute(arr);
        auto brute_force_time = watch.elapsed_time();

        // Print the results
        cout << left << setw(10) << size << setw(20) << optimized_time << setw(20) << brute_force_time << setw(20) << optimized_count << setw(20) << brute_force_count << endl;	
    }
}
// -----------------------------------------------------------------
// Function: run_tests_parallel
// This function measures the execution time of both the brute force
// and optimized "4-Sum" solutions in parallel for different input sizes.
//
// Version: 1.0
// -----------------------------------------------------------------
void run_tests_parallel(int start, int end, int step) {
    cout << left << setw(10) << "InputSize" << setw(20) << "OptimizedTime(ns)" << setw(20) << "BruteForceTime(ns)" << setw(20) << "OptimizedCount" << setw(20) << "BruteForceCount" << endl;	
    cout << string(90, '-') << endl;

    for (int size = start; size <= end; size += step) {
        vector<int> arr = create_list_of_numbers(size); // Generate input data

        StopWatch watch;

        // Future to get the result of optimized solution in a separate thread
        auto optimized_future = std::async(std::launch::async, [&arr, &watch]() {
            watch.reset(); // Start timer
            int optimized_count = four_sum(arr);
            return make_pair(optimized_count, watch.elapsed_time());
        });

        // Future to get the result of brute force solution in a separate thread
        auto brute_force_future = std::async(std::launch::async, [&arr, &watch]() {
            watch.reset(); // Start timer
            int brute_force_count = four_sum_brute(arr);
            return make_pair(brute_force_count, watch.elapsed_time());
        });

        // Wait for both computations to finish
        auto [optimized_count, optimized_time] = optimized_future.get();
        auto [brute_force_count, brute_force_time] = brute_force_future.get();

        // Print the results
        cout << left << setw(10) << size << setw(20) << optimized_time << setw(20) << brute_force_time << setw(20) << optimized_count << setw(20) << brute_force_count << endl;	
    }
}



void run_test_on_input_files(){
	std::map<int, std::string> fileMap = {
        {1000, "Assignment 1 Files/1Kints.txt"},
        {500, "Assignment 1 Files/5Hints.txt"},
        {8, "Assignment 1 Files/8ints.txt"},
        {50, "Assignment 1 Files/50ints.txt"},
        {100, "Assignment 1 Files/100ints.txt"},
        {200, "Assignment 1 Files/200ints.txt"}
    };
	cout << left << setw(10) << "InputSize" << setw(20) << "OptimizedTime(ns)" << setw(20) << "BruteForceTime(ns)" << setw(20) << "OptimizedCount" << setw(20) << "BruteForceCount" << endl;	

	//loop through the map
	for(auto const& x : fileMap){
		vector<int> arr = read_file(x.second);
		StopWatch watch;
		watch.reset();
		int optimized_count = four_sum(arr);
		auto optimized_time = watch.elapsed_time();
		watch.reset();
		int brute_force_count = four_sum_brute(arr);
		auto brute_force_time = watch.elapsed_time();
        cout << left << setw(10) << x.first << setw(20) << optimized_time << setw(20) << brute_force_time << setw(20) << optimized_count << setw(20) << brute_force_count << endl;	
	}

}


// -----------------------------------------------------------------
// Function: main
// The main function parses command-line arguments, then runs tests
// for the "4-Sum" problem across a range of input sizes.
// 
// Version: 1.0
// -----------------------------------------------------------------
int main(int argc, char* argv[]) {

	run_test_on_input_files();
    // if (argc != 4) { // Ensure correct usage
    //     cerr << "Usage: " << argv[0] << " <start> <end> <step>" << endl;
    //     return 1;
    // }

    // // Parse command-line arguments
    // int start = stoi(argv[1]);
    // int end = stoi(argv[2]);
    // int step = stoi(argv[3]);

    // run_tests_parallel(start, end, step); // Run tests with the specified range

    return 0;
}