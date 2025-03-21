// Huy Nguyen
// CS 6420
// Assignment #4
// Dr. Rague
// Due: 03/25/25
// Version: 1.0
// -----------------------------------------------------------------
// This program schedules activities using the Earliest Start Time algorithm.
// The activities are sorted based on the start time and scheduled on the first available resource.
// -----------------------------------------------------------------
// Usage: g++ -std=c++11 schedule.cpp -o schedule
// ./schedule
// or
// make schedule
// ./schedule
// -----------------------------------------------------------------


// Compiler directives
#include <iostream>
#include <vector>
#include <utility>
#include <tuple>
#include <random>
#include <iomanip>


// Function to check if two activities are compatible
// arguments: two activities
// returns: true if the activities are compatible, false otherwise
bool compatible(std::tuple<std::string, int, int> a, std::tuple<std::string, int, int> b){
	// compatible is when a is before b and a ends before b starts 
	if (std::get<2>(a) <= std::get<1>(b)){
		return true;
	}
	return false;
}


// Function to schedule activities using the Earliest Start Time algorithm
// arguments: vector of activities
// returns: vector of vectors of activities scheduled on each resourse

std::vector<std::vector<std::tuple<std::string, int, int>>> schedule_EST(std::vector<std::tuple<std::string, int, int>> activities){
	// This is Earliest Start Time algorithm to schedule activities
	// So sort the activities based on the start time
	// But assuming that the given activities are already sorted
	std::sort(activities.begin(), activities.end(), [](std::tuple<std::string, int, int> a, std::tuple<std::string, int, int> b){
		return std::get<1>(a) < std::get<1>(b);
	});

	// print the sorted activities
	// Find the longest element in the matrix to format the matrix printing
	std::cout << "Sorted activities: " << std::endl;
	int longest_element = 0;
	for(auto activity: activities){
		for (int i = 0; i < 3; i++){
			if (std::to_string(std::get<1>(activity)).length() > longest_element){
				longest_element = std::to_string(std::get<1>(activity)).length();
			}
		}
	}
	// print the activities
	for(auto activity: activities){
		//align left
		std::cout << std::setw(longest_element) << std::get<0>(activity) << " " << std::get<1>(activity) << " " << std::get<2>(activity) << std::endl;
	}
		
	// Schedule the activities
	// Create a vector of vectors to store the scheduled activities
	std::vector<std::vector<std::tuple<std::string, int, int>>> schedule_resourses;
	for(auto activity: activities){
		// check if there is a resourse available
		if (schedule_resourses.empty()){
			// if there is no resourse available, create a new resourse
			std::vector<std::tuple<std::string, int, int>> resourse;
			resourse.push_back(activity);
			schedule_resourses.push_back(resourse);
		}
		else{
			// check the current activity with all the resourses
			bool is_scheduled = false;
			for (auto it = schedule_resourses.rbegin(); it != schedule_resourses.rend(); ++it) {
				auto &resourse = *it;				
				// Because the activities are sorted, we only need to check the last activity in the resourse
				auto last_activity = resourse.back();
				if (compatible(last_activity, activity)){
					resourse.push_back(activity);
					is_scheduled = true;
					break;
				}
			}
			if (!is_scheduled){
				std::vector<std::tuple<std::string, int, int>> resourse;
				resourse.push_back(activity);
				schedule_resourses.push_back(resourse);
			}
		}
	}
	// print the schedule
	std::cout << "Optimal Schedule with " << schedule_resourses.size() << " resourses" << std::endl;
	for(int i=0; i<schedule_resourses.size(); i++){
		std::cout << "Resourse " << i+1 << ": ";
		for(auto activity: schedule_resourses[i]){
			std::cout << std::get<0>(activity) << " ";
		}
		std::cout << std::endl;
	}
	return schedule_resourses;
}


// Function to generate activities
// arguments: universal start time, universal end time, number of activities, maximum compatibility
// returns: vector of activities
std::vector<std::tuple<std::string, int, int>> generate_activities(int universal_start, int universal_end, int num_activities, int max_compatibility) {
    std::vector<std::tuple<std::string, int, int>> activities;
    std::vector<std::pair< int, int>> max_compatibility_array;

    std::random_device rd;
    std::mt19937 gen(rd());
	int chunk_size = (universal_end - universal_start) / num_activities;
	std::uniform_int_distribution<int> start_dist(universal_start, universal_end - chunk_size * max_compatibility/2);

	std::uniform_int_distribution<int> length_dist(chunk_size, chunk_size * max_compatibility/2);

	std::cout << "Chunk size: " << chunk_size << std::endl;
	int compatible_resourse_index = 0;
	bool new_resourse = true;
    for (int i = 1; i <= num_activities; i++) {
        int start, end;
        bool compatible = false;
		for(int j = 0; j < 20; j++) {
			start = start_dist(gen);
			end = start + length_dist(gen);
			compatible = false;
			
			for (int k = 0; k < max_compatibility_array.size(); k++) {
				auto activity = max_compatibility_array[k];
				int existing_start = activity.first;
				int existing_end = activity.second;

				bool is_compatible = (start >= existing_end || end <= existing_start);
				if(is_compatible) {
					compatible = true;
					compatible_resourse_index = k;
					break;
				}
			}
		}

		if(!compatible &&  max_compatibility_array.size() >= max_compatibility) {

			// pick a random activity from the max_compatibility_array
			std::uniform_int_distribution<int> activity_dist(0, max_compatibility_array.size() - 1);

			compatible_resourse_index = activity_dist(gen);
			
			
			start = max_compatibility_array.at(compatible_resourse_index).first;
			end = start + length_dist(gen);
			compatible = true;
			new_resourse = false;
		}

        std::string name = "a" + std::to_string(i);
        std::tuple<std::string, int, int> activity = {name, start, end};
        activities.push_back(activity);
		// add the activity to the max_compatibility_array with max_compatibility_array start and the new end
		if (max_compatibility_array.size() == 0) {
			max_compatibility_array.push_back({start, end});
		} else {
			if(new_resourse) {
				max_compatibility_array.push_back({start, end});
			} else {
				auto old_activity = max_compatibility_array[compatible_resourse_index];
				max_compatibility_array[compatible_resourse_index] = {old_activity.first, end};
			
			}
		}

        
    }
    
    return activities;
}

int main(int argc, char const *argv[]){
	
    // std::vector<std::tuple<std::string, int, int>> activities = {
    //     {"a1", 900, 1030}, {"a2", 900, 1230}, {"a3", 900, 1030}, 
    //     {"a4", 1100, 1230}, {"a5", 1100, 1400}, {"a6", 1300, 1430}, 
    //     {"a7", 1300, 1430}, {"a8", 1400, 1630}, {"a9", 1500, 1630}, {"a10", 1500, 1630}
	// };
	// set 1:
	// universal_start = 200
	// universal_end = 1200
	// num_activities = 6
	// max_compatibility = 2

	// set 2:
	// universal_start = 1000
	// universal_end = 2000
	// num_activities = 10
	// max_compatibility = 3

	// set 3:
	// universal_start = 800
	// universal_end = 1500	
	// num_activities = 20
	// max_compatibility = 4

	// set 4:
	// universal_start = 500
	// universal_end = 2000
	// num_activities = 30
	// max_compatibility = 5

	// set 5:
	// universal_start = 700
	// universal_end = 1800
	// num_activities = 50
	// max_compatibility = 7

	std::vector<std::tuple<int, int, int, int>> sets = {
		{900, 1200, 6, 3},
		{1000, 2000, 10, 3},
		{800, 1500, 20, 4},
		{500, 2000, 30, 5},
		{700, 1800, 50, 7}
	};

	for (int i = 0; i < sets.size(); i++) {
		std::vector<std::tuple<std::string, int, int>> activities = generate_activities(std::get<0>(sets[i]), std::get<1>(sets[i]), std::get<2>(sets[i]), std::get<3>(sets[i]));
		std::cout << "Test: " << i + 1 << std::endl << std::endl;
		schedule_EST(activities);
		std::cout << std::endl;
	}

	return 0;
}