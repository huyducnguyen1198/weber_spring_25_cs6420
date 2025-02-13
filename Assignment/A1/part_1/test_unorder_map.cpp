#include <iostream>
#include <unordered_map>
#include <array>
#include <chrono>
#include <iomanip>

struct ArrayHash {
    std::size_t operator()(const std::array<int, 4>& arr) const {
        std::hash<int> hasher; // Create a hash function object
        std::size_t h1 = hasher(arr[0]);
        std::size_t h2 = hasher(arr[1]);
        std::size_t h3 = hasher(arr[2]);
        std::size_t h4 = hasher(arr[3]);
        return h1 ^ (h2 << 1) ^ (h3 << 2) ^ (h4 << 3);
    }
};

void testUnorderedMapComplexity(int start, int end, int step) {
    std::cout << std::setw(10) << "N" << std::setw(20) << "Time (ms)" << "\n";
    std::cout << std::string(30, '-') << "\n";

    for (int n = start; n <= end; n += step) {
        std::unordered_map<std::array<int, 4>, int, ArrayHash> quad_map;

        auto start_time = std::chrono::high_resolution_clock::now();

        for (int i = 0; i < n; ++i) {
            std::array<int, 4> temp = {i, i * 2, i * 3, i * 4};
            if (quad_map.find(temp) == quad_map.end()) {
                quad_map[temp] = 1;
            }
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();

        std::cout << std::setw(10) << n << std::setw(20) << duration << "\n";
    }
}

int main() {
    int start, end, step;

    std::cout << "Enter the start value for N: ";
    std::cin >> start;

    std::cout << "Enter the end value for N: ";
    std::cin >> end;

    std::cout << "Enter the step size: ";
    std::cin >> step;

    testUnorderedMapComplexity(start, end, step);

    return 0;
}
