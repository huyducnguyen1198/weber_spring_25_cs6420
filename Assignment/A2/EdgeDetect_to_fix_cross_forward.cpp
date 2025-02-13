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
// To run: ./EdgeDetect graph.txt
// To clean: make clean

// Compiler directives
#include "alg_graphs.h"
#include <fstream>



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

	DepthFirstSearch dfs(graph);

	dfs.print_edge_types();


}