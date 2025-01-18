#include <iostream>
#include <vector>
#include <unordered_map>	
#include <algorithm>
#include <set>
#include <fstream>
#include <sstream>
#include "alg_stopwatch.h"
#include <unordered_set>
#include <map>
#include <array>
#include <list>
using namespace std;

int four_sum_brute(const vector<int> & arr){
	int total_pairs = 0;
	set<vector<int> > pairs;

	for(int i = 0; i < arr.size(); i++){
		for(int j = i+1; j < arr.size(); j++){
			for(int k = j+1; k < arr.size(); k++){
				for(int l = k+1; l < arr.size(); l++){
					if(arr[i] + arr[j] + arr[k] + arr[l] == 0){
						vector<int> temp = {arr[i], arr[j], arr[k], arr[l]};
						sort(temp.begin(), temp.end());
						pairs.insert(temp);
						total_pairs++;
					}
				}
			}
		}
	}

	return total_pairs;
}

void sort_four(int& a, int& b, int& c, int& d) {
    if (a > b) std::swap(a, b);
    if (c > d) std::swap(c, d);
    if (a > c) std::swap(a, c);
    if (b > d) std::swap(b, d);
    if (b > c) std::swap(b, c);
}


struct ArrayHash {
    template <typename T, std::size_t N>
    std::size_t operator()(const std::array<T, N>& arr) const {
        std::size_t hash = 0;
        for (const auto& elem : arr) {
            hash ^= std::hash<T>{}(elem) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};


int four_sum(const vector<int> & arr){

	// Get all the pairs, and store them in a map with the sum as the key
	unordered_map<int, vector<pair<int, int> > > sum_pairs;
	int total_pairs = 0;
	for(int i = 0; i < arr.size(); i++){
		for(int j = i+1; j < arr.size(); j++){
			sum_pairs[arr[i] + arr[j]].push_back({arr[i], arr[j]});
			total_pairs++;
		}
	}
	
	//linked list
	
    std::unordered_map<std::array<int, 4>, std::string, ArrayHash> quad_map;

	std::unordered_map<std::string, int> quad_map_str;
	quad_map.max_load_factor(.6);
	quad_map.rehash(total_pairs);
	quad_map.reserve(total_pairs);

	for(const auto & [sum, pairs]: sum_pairs){
		if(sum_pairs.find(-sum) != sum_pairs.end()){
			for(const auto & [a, b]: pairs){
				for(const auto & [c, d]: sum_pairs[-sum]){
					if(a != c && a != d && b != c && b != d){
						// Updated Code Using sort_four
						int a_sorted = a, b_sorted = b, c_sorted = c, d_sorted = d;
						sort_four(a_sorted, b_sorted, c_sorted, d_sorted);

						// std:string temp_str = to_string(a_sorted) + "_" + to_string(b_sorted) + "_" + to_string(c_sorted) + "_" + to_string(d_sorted);
						// if (quad_map_str.find(temp_str) == quad_map_str.end()) {
						// 	quad_map_str[temp_str] = 1;
						// }

						array<int, 4> temp = {a_sorted, b_sorted, c_sorted, d_sorted};
						if (quad_map.find(temp) == quad_map.end()) {
							quad_map[temp] = 1;
						}
						// create vector on heap
					}
				}
			}
			
		}

	}
	return quad_map.size();
}



vector<int> read_file(const string& filename) {
    /**
     * Function to read file and return a vector of integers
     * If the first line contains more than one value, read all values from the first line.
     * If the first line contains only one value, read each subsequent value from separate lines.
     * @param filename: The name of the file to read
     * @return: A vector of integers
     */

    ifstream file(filename);
    vector<int> arr;

    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return arr;
    }

    string line;
    // Read the first line
    if (getline(file, line)) {
        istringstream iss(line);
        int num;
        vector<int> first_line_values;

        while (iss >> num) {
            first_line_values.push_back(num);
        }

        // If more than one value is found in the first line, add them all and stop reading further
        if (first_line_values.size() > 1) {
            arr.insert(arr.end(), first_line_values.begin(), first_line_values.end());
        } else if (!first_line_values.empty()) {
            // If only one value is found, continue reading one value per line
            arr.push_back(first_line_values[0]);

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
vector<int> read_file_multiple_lines(const string& filename) {
    /**
     * Function to read file and return a vector of integers
     * The file is assumed to have integers spread across multiple lines
     * @param filename: The name of the file to read
     * @return: A vector of integers
     */

    ifstream file(filename);
    vector<int> arr;

    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return arr;
    }

    string line;
    while (getline(file, line)) {
        istringstream iss(line);
        int num;
        while (iss >> num) {
            arr.push_back(num);
        }
    }

    file.close();
    return arr;
}
vector<int> generate_random_unique_array(int size) {
    /**
     * Generate an array of random unique integers.
     * @param size: The number of unique integers to generate
     * @return: A vector of random unique integers
     */
    vector<int> arr;
    unordered_set<int> unique_values;

    srand(static_cast<unsigned>(time(0))); // Seed for random number generation

    while (unique_values.size() < size) {
        int random_value = rand();
        unique_values.insert(random_value);
    }

    arr.assign(unique_values.begin(), unique_values.end());
    return arr;
}

vector<int> create_list_of_numbers(int n) {
    /**
     * Generate a vector of random integers in the range [-n, n].
     * @param n: The range of numbers and the size of the list to generate
     * @return: A vector of random integers
     */
    vector<int> arr;
    srand(static_cast<unsigned>(time(0))); // Seed for random number generation

    for (int i = 0; i < n; ++i) {
        int random_value = rand() % (2 * n + 1) - n; // Generate random value in range [-n, n]
		if (find(arr.begin(), arr.end(), random_value) != arr.end()) {
			// If the value already exists in the array, generate a new value
			--i;
			continue;
		}
        arr.push_back(random_value);
    }

    return arr;
}
// Function to run tests and measure execution times
void run_tests(int start, int end, int step) {
    cout << left << setw(10) << "InputSize" << setw(20) << "OptimizedTime(ns)" << setw(20) << "BruteForceTime(ns)" << endl;// setw(20) << "OptimizedCount" << setw(20) << "BruteForceCount" << endl;
    cout << string(100, '-') << endl;
	StopWatch watch;

    for (int size = start; size <= end; size += step) {
        // Generate an array of random unique values
        vector<int> arr = create_list_of_numbers(size);
		


        // Measure optimized function time
        watch.reset();
        int optimized_count = four_sum(arr);
        auto optimized_time = watch.elapsed_time();

        watch.reset();
        int brute_force_count = four_sum_brute(arr);
        auto brute_force_time = watch.elapsed_time();

        cout << left << setw(10) << size << setw(20) << optimized_time << setw(20) << brute_force_time<<endl;// << setw(20) << optimized_count << setw(20) << brute_force_count << endl;

	}
}
void run_test_on_input_files(){
	std::map<int, std::string> fileMap = {
        {1000, "Assignment/A1/part_1/Assignment 1 Files/1Kints.txt"},
        {500, "Assignment/A1/part_1/Assignment 1 Files/5Hints.txt"},
        {8, "Assignment/A1/part_1/Assignment 1 Files/8ints.txt"},
        {50, "Assignment/A1/part_1/Assignment 1 Files/50ints.txt"},
        {100, "Assignment/A1/part_1/Assignment 1 Files/100ints.txt"},
        {200, "Assignment/A1/part_1/Assignment 1 Files/200ints.txt"}
    };

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
		cout << left << setw(10) << x.first << setw(20) << optimized_time << setw(20) << brute_force_time << 
	}

}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        cerr << "Usage: " << argv[0] << " <start> <end> <step>" << endl;
        return 1;
    }

    // Parse command-line arguments
    int start = stoi(argv[1]);
    int end = stoi(argv[2]);
    int step = stoi(argv[3]);

    // Run the tests
    run_tests(start, end, step);

    return 0;
}
