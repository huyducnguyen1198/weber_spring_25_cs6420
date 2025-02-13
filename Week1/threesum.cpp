#include <iostream>
#include <vector>
#include <unordered_map>
#include <set>
#include <cstdlib>
#include <ctime>
#include <algorithm>

// Function to create a list of random numbers
std::vector<int> create_list_of_numbers(int n) {
    std::vector<int> arr;
    for (int i = 0; i < n; ++i) {
        arr.push_back(rand() % 11001 - 1000); // Random number between -1000 and 10000
    }
    return arr;
}

// Optimized three_sum function
int three_sum(const std::vector<int>& arr) {
    std::unordered_map<int, std::vector<std::pair<int, int>>> sum_of_two;
    for (size_t i = 0; i < arr.size(); ++i) {
        for (size_t j = i + 1; j < arr.size(); ++j) {
            int suum = arr[i] + arr[j];
            sum_of_two[suum].emplace_back(arr[i], arr[j]);
        }
    }

    std::set<std::tuple<int, int, int>> sum_of_three;
    for (size_t i = 0; i < arr.size(); ++i) {
        int target = -arr[i];
        if (sum_of_two.find(target) != sum_of_two.end()) {
            for (const auto& pair : sum_of_two[target]) {
                if (arr[i] != pair.first && arr[i] != pair.second) {
                    std::vector<int> triplet = {arr[i], pair.first, pair.second};
                    std::sort(triplet.begin(), triplet.end());
                    sum_of_three.emplace(triplet[0], triplet[1], triplet[2]);
                }
            }
        }
    }

    return sum_of_three.size();
}

// Brute force three_sum function
int three_sum_brute(const std::vector<int>& arr) {
    int total_pairs = 0;
    std::set<std::vector<int>> all_pairs;
    for (size_t i = 0; i < arr.size(); ++i) {
        for (size_t j = i + 1; j < arr.size(); ++j) {
            for (size_t k = j + 1; k < arr.size(); ++k) {
                if (arr[i] + arr[j] + arr[k] == 0) {
                    std::vector<int> triplet = {arr[i], arr[j], arr[k]};
                    std::sort(triplet.begin(), triplet.end());
                    if (all_pairs.find(triplet) == all_pairs.end()) {
                        all_pairs.insert(triplet);
                        ++total_pairs;
                    }
                }
            }
        }
    }

    return total_pairs;
}

// Test and print results
void test_and_print_three_sum() {
    std::srand(std::time(0)); // Seed random number generator
    std::vector<int> sizes;
    for (int n = 10; n <= 5000; n += 200) {
        sizes.push_back(n);
    }

    std::vector<double> optimized_times;
    std::vector<double> brute_force_times;

    for (int size : sizes) {
        std::vector<int> arr = create_list_of_numbers(size);

        // Measure optimized time
        clock_t start_time = clock();
        int num_pair_opt = three_sum(arr);
        clock_t end_time = clock();
        optimized_times.push_back(double(end_time - start_time) / CLOCKS_PER_SEC);

        // Measure brute force time (for small sizes only to avoid timeout)
        start_time = clock();
        int num_pair_brute =  three_sum_brute(arr);
        end_time = clock();
        brute_force_times.push_back(double(end_time - start_time) / CLOCKS_PER_SEC);

        std::cout << "Size: " << size
                  << " | Optimized Time: " << optimized_times.back()
                  << "s | Brute Force Time: " << brute_force_times.back()
				  << "s | Optimized Pairs: " << num_pair_opt
				  << " | Brute Force Pairs: " << num_pair_brute
				  << std::endl;
    }
}

int main() {
    test_and_print_three_sum();
    return 0;
}
