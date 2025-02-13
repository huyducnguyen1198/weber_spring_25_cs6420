// Huy Nguyen
// CS 6420
// Assignment #2
// Dr. Rague
// Due: 01/31/25
// Version: 1.0
// -----------------------------------------------------------------
// This program reads a directed graph from a file and performs
// topological sorting using Kahn’s algorithm, which is a BFS. It also checks if
// the topological order is unique.
// -----------------------------------------------------------------
// Part 1:
// The time complexity of the Kahn’s algorithm is O(V + E) where V is the number of vertices
// and E is the number of edges. The algorithm iterates through all the vertices and edges.
// The loop eventually add all the vertices with 0 indegree to the zero_indegree set, 
// it does so by iterating and removing all the edges. 
// The act of inserting and removing from the set is O(logV) because of the set data structure(red-black tree).
// So the time complexity of the entire loop is O(VlogV + E).
// The function to check the unique topological order is the same as the original Kahn's algorithm 
// except for the addition of the vertex_color map, which is O(1) to insert and access.

// Part 2:
// The provided code is a recursive DFS, which has the time complexity of O(V + E) where V is the number of vertices
// and E is the number of edges. 
// Additions are made under '/ /added' comments in the provided code.
// THe core structure is the addition of vector of tuples to store the edge types and the function to print them.
// Vector has O(1) complexity to insert and access.
// So the time complexity of the entire loop is O(V + E).


// To complie: make
// To run: ./BFSTopo graph.txt
// To clean: make clean








// Compiler directives
#include "alg_graphs.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <unordered_map>

// -----------------------------------------------------------------
// This function performs Kahn’s algorithm for topological sorting.
//
// Parameters:
// - graph: A directed graph
//
// If a cycle is detected, an error message is printed, and the
// program exits.
// DONT NEED THIS FUNCTION. IT IS JUST FOR REFERENCE
// -----------------------------------------------------------------
void kahn_algo(Digraph graph){
	
	int num_vertices = graph.BaseGraph::V(); // O(1)
	std::set <int> zero_indegree; // O(1)

	std::vector<int> topological_order; // O(1)

	for (int i = 0; i < num_vertices; i++){ // O(VlogV) Set insert is O(logV)
		if(graph.in_degree(i) == 0){
			zero_indegree.insert(i);
		}
	}

	while (!zero_indegree.empty()){ // O(V + E) but set insert and erase is O(logV) => O(VlogV + E)
		// Check one of the vertex with 0 indegree
		int vertex = *zero_indegree.begin();
		zero_indegree.erase(zero_indegree.begin());

		// Add the vertex to the topological order
		topological_order.push_back(vertex);

		for (int v : graph.adj(vertex)){
			graph.remove_edge(vertex, v);

			if(graph.in_degree(v) == 0){
				zero_indegree.insert(v);
			}
		}
	}

	
	if(topological_order.size() != num_vertices){
		std::cerr << "Error: Graph has a cycle" << std::endl;
		exit(1);
	}

	std::cout << "Topological order: ";
	for (int i = 0; i < topological_order.size(); i++){
		std::cout << topological_order[i] << " ";
	}
}


// -----------------------------------------------------------------
// This function performs Kahn’s algorithm and determines if the
// topological order is unique.
//
// Parameters:
// - graph: A directed graph
//
// The function modifies Kahn’s algorithm by tracking visited nodes
// to determine uniqueness.
// -----------------------------------------------------------------

void kahn_algo_update_for_unique(Digraph graph){
	
	// Initialize variables
	int num_vertices = graph.BaseGraph::V(); 				// O(1)
	
	// Initialize data structures
	// Zero indegree vertices used to track vertices with 0 indegree
	
	std::set <int> zero_indegree; 							// O(1)

	// Topological order used to store the topological order
	std::vector<int> topological_order;						// O(1)

	// Vertex color used to track visited vertices
	// 0: white, 1: grey, 2: black
	std::unordered_map<int, int> vertex_color; 				// O(1)		




	// Initialize zero_indegree set and vertex_color map
	// by iterating through all vertices 
	// and adding vertices with 0 indegree to zero_indegree set
	for (int i = 0; i < num_vertices; i++){ 				// O(VlogV) Set insert is O(logV)
		if(graph.in_degree(i) == 0){
			zero_indegree.insert(i);         	  
		}
		vertex_color[i] = 0; 								// O(1) unordered_map insert is O(1) because of hashing
	}

	// for (int i : zero_indegree){
	// 		std::cout << i << " ";
	// }
	// std::cout << std::endl;


	// Initialize is_unique to true
	bool is_unique = true;									// O(1)		

	// Loop through the zero_indegree set
	// Which will have all the vertices 
	// Because each edge met will be removed, 
	// and the vertex will eventually have 0 indegree 
	// This will be O(V) because we are iterating through all the vertices
	// in conclusion, the time complexity of this entire loop is O(V + E)
	// but because of the usage of Set, and each vertex is only added once and removed once 
	// so the time complexity is O(VlogV + E)
	while (!zero_indegree.empty()){							
		// Check one of the vertex with 0 indegree 
		int vertex = *zero_indegree.begin(); 				// O(1)

		// Remove the vertex from the zero_indegree set
		zero_indegree.erase(zero_indegree.begin());			// O(logV) Set erase is O(logn) potentiall n = V
		

		// check color of vertex
		// if the vertex is visited_vertex and still has an indegree
		// meaning it can be reached by another vertex
		// the graph has multiple topological orders, set is_unique to false
		if(vertex_color[vertex] != 0 && graph.in_degree(vertex) > 0){ // O(1)
			is_unique = false;
		}

		
		// Add the vertex to the topological order
		// because it has 0 indegree
		topological_order.push_back(vertex);				// O(1)

		// Change the color of the vertex to grey
		vertex_color[vertex] = 1; 							// O(1)


		// Loop through the adjacency list of the vertex
		// to remove the edges and update the indegree of the vertices
		// This will be O(E) because we are iterating through all the edges
		for (int v : graph.adj(vertex)){  

			// check if v is visited_vertex
			if(vertex_color[v] == 0){								// O(1)	
				// Change the color of the vertex to grey
				vertex_color[v] = 1;										
			}else if(vertex_color[v] != 0 && graph.in_degree(v) > 0){  // O(1)
				// Again if it is visited_vertex and still has an indegree
				// meaning it can be reached by another vertex
				// This is only reached when the vertex is visited_vertex
				// but it might have more than 1 indegree
				is_unique = false;
			}

			// Remove the edge from the graph
			graph.remove_edge(vertex, v); 					// O(1) because of linked list

			// Check if the vertex has 0 indegree
			if(graph.in_degree(v) == 0){
				zero_indegree.insert(v);					// O(logV) Set insert is O(logV)
			}
			

			
			// std::cout << "checking from " << vertex << " to " << v << std::endl;

			// std::cout << graph<< std::endl;

			// std::cout << "Is unique: " << is_unique << std::endl;
		}

		// Change the color of the vertex to black
		vertex_color[vertex] = 2; 							// O(1)
	}

	// Check if the graph has a edge, it has a cycle
	if(graph.BaseGraph::E() > 0){ 							
		std::cerr << "Error: Graph has a cycle" << std::endl;
		exit(1);
	}
	
	// Print the topological order
	std::cout << "Topological sort: "; 
		for (int i = 0; i < topological_order.size(); i++){
			std::cout << topological_order[i] << " ";
		}
	std::cout << std::endl;

	// Check if the graph has multiple topological orders
	std::cout << "Unique: ";
	if(is_unique){

		std::cout << "Yes" << std::endl;
	}else{
		
		std::cout << "No" << std::endl;
	}
}



// -----------------------------------------------------------------
// Main function that reads the graph and runs Kahn's algorithm.
//
// Parameters:
// - argc: Number of command-line arguments
// - argv: Array of command-line arguments
//
// The program expects a single argument specifying the graph file path.
// -----------------------------------------------------------------

int main(int argc, char *argv[]){

	// Check if the correct number of arguments are provided
	if (argc != 2){
		std::cerr << "Usage: " << argv[0] << " <graph file>" << std::endl;
		exit(1);
	}

	// Get the graph file path from the command line arguments
	std::string graph_path = argv[1];

	// Open the graph file
	std::ifstream file(graph_path);

	// Declare a Digraph object
	Digraph graph;

	// Read the graph from the file and store it in the Digraph object
	file >> graph;

	// Perform Kahn’s algorithm for topological sorting
	// This function will print the topological order
	// and exit if a cycle is detected
	kahn_algo_update_for_unique(graph);

}