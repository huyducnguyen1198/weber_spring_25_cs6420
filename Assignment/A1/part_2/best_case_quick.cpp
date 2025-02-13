// Huy Nguyen
// CS 6420
// Assignment #1
// Dr. Rague
// Due: 01/16/25
// Version: 1.0
// -----------------------------------------------------------------
// This program generates a best-case input for quicksort using recursion and iteration.
// It includes testing functionalities to validate the output for both methods.
// Usage: ./best_case_quick 10
// Compilation: g++ -std=c++17 best_case_quick.cpp -o best_case_quick
// -----------------------------------------------------------------

#include <iostream>   // For input and output operations
#include <vector>     // For using the std::vector container
#include <iomanip>    // For formatted output (e.g., std::setw)
#include <string>     // For std::string
#include <utility>    // For std::pair
#include <stack>      // For std::stack

// -------------------------------------------------------------------
// Generates a sorted input string for testing quicksort.
// Example: n=4 produces "ABCD".
// -------------------------------------------------------------------
std::string generate_sorted_input_string(int n) {
    std::string input_str = ""; // Initialize an empty string
    for (int i = 0; i < n; i++) {
        input_str += 'A' + i; // Append the i-th letter of the alphabet
    }
    return input_str; // Return the generated string
}

// -------------------------------------------------------------------
// Recursively partitions an array for the best-case quicksort.
// Returns the partitioned array.
// -------------------------------------------------------------------
std::vector<char> partition(std::vector<char> arr) {
    int start = 0; // Start index
    int end = arr.size(); // End index

    // Base case: if one or no elements, return the array
    if (end - start <= 1) return arr;

    int mid = end % 2 == 0 ? end / 2 - 1 : end / 2; // Calculate the mid index

    // Create pivot, left, and right partitions
    std::vector<char> pivot(1, arr[mid]);
    std::vector<char> left(arr.begin() + start, arr.begin() + mid);
    std::vector<char> right(arr.begin() + mid + 1, arr.end());

    // Recursively partition the left and right subarrays
    std::vector<char> combined = partition(left); // Partition the left subarray
    combined.insert(combined.end(), pivot.begin(), pivot.end()); // Add the pivot
    std::vector<char> rightPartition = partition(right); // Partition the right subarray
    combined.insert(combined.end(), rightPartition.begin(), rightPartition.end()); // Add the right subarray

    // Swap the pivot with the start element
    std::swap(combined[start], combined[start + (mid - start)]);
    return combined; // Return the combined array
}

// -------------------------------------------------------------------
// In-place recursive partitioning of an array for quicksort.
// Operates using indices instead of subarrays.
// -------------------------------------------------------------------
void partition_v2(std::vector<char>& arr, int start, int end) {
    if (end - start <= 1) return; // Base case: single element

    int mid = ((end - start) % 2 != 0) ? ((end - start) / 2) : ((end - start) / 2 - 1); // Mid index
    mid = start + mid; // Adjust mid relative to the start

    // Recursively partition left and right subarrays
    partition_v2(arr, start, mid);
    partition_v2(arr, mid + 1, end);

    // Swap the pivot with the start element
    std::swap(arr[mid], arr[start]);
}

// -------------------------------------------------------------------
// Generates the best-case order for quicksort using recursion.
// -------------------------------------------------------------------
std::string best_case_quick_sort_rec(std::string input_str) {
    std::vector<char> arr(input_str.begin(), input_str.end()); // Convert string to vector
    partition_v2(arr, 0, arr.size()); // Partition the array
    return std::string(arr.begin(), arr.end()); // Convert back to string
}

// -------------------------------------------------------------------
// Generates the best-case order for quicksort using iteration.
// -------------------------------------------------------------------
std::string best_case_quick_sort_iter(std::string input_str) {
    std::vector<char> arr(input_str.begin(), input_str.end()); // Convert string to vector
    std::stack<std::pair<int, int>> pair_stack; // Stack for iterative partitioning
    std::vector<std::pair<int, int>> to_swap; // Track elements to swap

    // Push initial indices to the stack
    pair_stack.push({0, arr.size()});

    // Perform iterative partitioning
    while (!pair_stack.empty()) {
        std::pair<int, int> current_pair = pair_stack.top(); // Get top indices
        int start = current_pair.first;
        int end = current_pair.second;
        pair_stack.pop(); // Remove top of stack

        int mid = ((end - start) % 2 != 0) ? ((end - start) / 2) : ((end - start) / 2 - 1); // Mid index
        mid += start; // Adjust mid relative to start

        // Push left and right partitions if valid
        if (end - start > 1) {
            pair_stack.push({start, mid});
            pair_stack.push({mid + 1, end});
            to_swap.push_back({start, mid}); // Record the swap
        }
    }

    // Perform swaps in reverse order
    for (int i = to_swap.size() - 1; i >= 0; i--) {
        std::pair<int, int> current_pair = to_swap[i];
        int start = current_pair.first;
        int end = current_pair.second;
        std::swap(arr[start], arr[end]); // Swap elements
    }

    return std::string(arr.begin(), arr.end()); // Convert back to string
}

// -------------------------------------------------------------------
// Utility function to count the number of digits in an integer.
// -------------------------------------------------------------------
int int_letter(int n) {
    return std::to_string(n).size(); // Convert to string and return length
}

// -------------------------------------------------------------------
// Main function to test and compare recursive and iterative methods.
// -------------------------------------------------------------------
int main(int argv, char* argc[]) {
    // Validate command-line arguments
    if (argv != 2) {
        std::cout << "Usage: " << argc[0] << " <n>\n";
        return 1;
    }

    int n = std::stoi(argc[1]); // Parse input size

    std::vector<std::string> rec_outputs; // Store recursive outputs
    std::vector<std::string> iter_outputs; // Store iterative outputs

    // Print table header
    int num_letter_n = int_letter(n);
    std::cout << "n" << std::setw(n + 2) << "input" << std::setw(n + 1) << "rec_output"
              << std::setw(n + 1) << "iter_output" << std::setw(6) << "equal" << std::endl;

    // Loop through sizes and test outputs
    for (int i = 1; i < n; i++) {
        std::string input_str = generate_sorted_input_string(i); // Generate input
        std::string output_rec = best_case_quick_sort_rec(input_str); // Recursive output
        std::string output_iter = best_case_quick_sort_iter(input_str); // Iterative output
        bool is_equal = output_rec == output_iter; // Check equality

        int num_letter_i = int_letter(i); // Count digits in `i`

        // Print row of results
        std::cout << i << std::setw(n + 2 - num_letter_i) << input_str
                  << std::setw(n + 1) << output_rec << std::setw(n + 1) << output_iter
                  << std::setw(6) << is_equal << std::endl;
    }

    return 0; // Exit program
}