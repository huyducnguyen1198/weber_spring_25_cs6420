/******************************************************************************
 *  File: alg_graphs.cpp
 * 
 *  An implementation file of multiple graph-related classes. 
 *
 ******************************************************************************/

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include "alg_graphs.h"

/******************************************************************************
 *  Class: BaseGraph
 *  A base graph class capturing what is common between undirected and directed
 *  classes.
 ******************************************************************************/
//checks that a vertex is contained in the graph
void BaseGraph::validate_vertex(int v) const { 
  if(v < 0 || v > _V) {
    throw std::runtime_error("vertex " + std::to_string(v) + " is not between 0 and " + std::to_string(_V - 1));
  }
}

//copies a graph onto another graph
void BaseGraph::copy_graph(const BaseGraph &g) {
  V(g.V());
  _E = g.E();

  for (int v = 0; v < _V; v++) {
    for(int w : g.adj(v)) {
      _adj[v].push_back(w);
    }
  }
}

// Constructors
// Creates a graph, and sets the adjacency list to contain V vertices
BaseGraph::BaseGraph(int V): _V(V) {
  if(_V != 0){
    _adj = new std::list<int>[V];
  }
}

// copy constructor and operator, uses the copy method
BaseGraph::BaseGraph(const BaseGraph& g): _V(0){
  copy_graph(g);
}

// = overloaded operator for a graph
BaseGraph& BaseGraph::operator=(const BaseGraph& g){
  if(this != &g) {
    delete[] _adj;
    _V = 0;
    copy_graph(g);
  }
  
  return *this;
}

// Move constructor and operator
// constructs graph using another graph, then resets the contents of the other graph
BaseGraph::BaseGraph(BaseGraph&& g) noexcept: _V(g._V), _E(g._E), _adj(g._adj) {
  g._V = 0;
  g._E = 0;
  g._adj = nullptr;
}

// another overloaded = operator
BaseGraph& BaseGraph::operator=(BaseGraph&& g) noexcept{
  delete[] _adj;
  _V = g._V;
  _E = g._E = 0;
  _adj = g._adj;

  g._V = 0;
  g._E = 0;
  g._adj = nullptr;

  return *this;
}

// Vertices and edges

// returns the number of vertices in a base graph
int BaseGraph::V() const { return _V; }


// creates a new adjacency list for a not yet made graph?
void BaseGraph::V(int V) { 
  if(_V != 0 && _V != V) {
    throw std::runtime_error("Cannot change the number of vertices of an existing graph");
  }

  _V = V; 
  _adj = new std::list<int>[V];
}

// returns number of edges
int BaseGraph::E() const { return _E; }

// returns true if an edge exists between two vertices
bool BaseGraph::edge(int v, int w) const {
  validate_vertex(v);
  validate_vertex(w);
  auto &list = _adj[v];
  auto it = std::find(list.begin(), list.end(), w);
  return it != list.end();
}

// returns the adjacency list of vertex v
std::list<int> BaseGraph::adj(int v) const {
  validate_vertex(v);
  return _adj[v];
}

// Input/output
std::string BaseGraph::str() const {
  std::ostringstream sout;
  sout << *this;
  return sout.str(); 
}

std::ostream& operator<<(std::ostream &out, const BaseGraph& g){
  out << g._V << std::endl << g._E << std::endl;
  for (int v = 0; v < g._V; v++) {
    out << v << ": ";
    for (int w : g._adj[v]) {
      out << w << " ";
    }

    out << std::endl;
  }

  return out; 
}

std::istream& operator>>(std::istream &in, BaseGraph& g){
  std::string line;
  while(std::getline(in,line)){
    if(g._V == 0){
      g.V(std::stoi(line));
    } else {
      std::istringstream ss(line);

      int v; ss >> v;
      char colon; ss >> colon;
      int w;
      while(ss >> w){
        if(g.is_directed()){
          g.add_edge(v, w);
        } else {
          g._E++;
          g._adj[v].push_back(w);
        }
      }
    }
  }

  if(!g.is_directed()) g._E /= 2;

  return in;
}

// Clean up
BaseGraph::~BaseGraph() noexcept{
  delete[] _adj;
}

/******************************************************************************
 *  Class: Graph
 *  A class representing undirected graphs
 ******************************************************************************/
// Constructors
Graph::Graph(int V): BaseGraph(V){}

  // Copy and move constructors
Graph::Graph(Graph& g): BaseGraph(g){}
Graph::Graph(Graph&& g) noexcept: BaseGraph(g){}

// returns false saying that the graph is undirected
bool Graph::is_directed() const { return false; }

// returns the degree of a vertex (returns the adjacency list size)
int Graph::degree(int v) const {
  return this->adj(v).size();
}

// Adding/removing edges

// checks that two vertices exist, and if so, creates an edge (adds each other to their adjacency lists)
void Graph::add_edge(int v, int w) {
  validate_vertex(v);
  validate_vertex(w);
  _E++;
  _adj[v].push_back(w);
  _adj[w].push_back(v);
}

// checks that two vertices exist, and if there is an edge between them, it removes it
void Graph::remove_edge(int v, int w) {
  validate_vertex(v);
  validate_vertex(w);
  auto &v_list = this->_adj[v];
  auto it = std::find(v_list.begin(), v_list.end(), w);
  v_list.erase(it);

  auto w_list = this->_adj[w];
  it = std::find(w_list.begin(), w_list.end(), v);
  w_list.erase(it);
  this->_E--;
}

// Cleanup
Graph::~Graph() noexcept {}

/******************************************************************************
 *  Class: Digraph
 *  A class representing directed graphs
 ******************************************************************************/
// Constructors
Digraph::Digraph(int V): BaseGraph(V) {
  if(V != 0){
    indegree = new int[V]{0}; 
  }
}

// Copy constructor and  assignment operator
Digraph::Digraph(const Digraph& g): BaseGraph(g){
  if(this->_V != 0){
    indegree = new int[this->_V]{0}; 
    for(int v = 0; v < this->_V; v++) {
      indegree[v] = g.indegree[v];
    }
  }
}

Digraph& Digraph::operator=(const Digraph& g){
  if(this != &g) {
    this->BaseGraph::operator=(g);

    if(this->_V != 0){
      indegree = new int[this->_V]{0}; 
      for(int v = 0; v < this->_V; v++) {
        indegree[v] = g.indegree[v];
      }
    }
  }
  
  return *this;
}

// Move constructor and  assignment operator
Digraph::Digraph(Digraph&& g) noexcept: BaseGraph(g) {
  indegree = g.indegree;
  g.indegree = nullptr;
}

Digraph& Digraph::operator=(Digraph&& g) noexcept {
  this->BaseGraph::operator=(g);
  indegree = g.indegree;
  g.indegree = nullptr;

  return *this;
}

// Vertices and edges

// creates a digraph with V vertices, and sets the in degree of all vertices to 0
void Digraph::V(int V) { 
  BaseGraph::V(V);
  indegree = new int[V]{0}; 
}

// returns true saying the graph is directed
bool Digraph::is_directed() const { return true; }

// Degrees

// returns the degree of a vertex, adds in_degree and out_degree
int Digraph::degree(int v) const {
  return out_degree(v) + in_degree(v);
}

// returns the out degree of a vertex (returns the size of its adjacency list)
int Digraph::out_degree(int v) const {
  return this->adj(v).size();
}

// returns the in degree of a vertex
int Digraph::in_degree(int v) const {
  this->validate_vertex(v);
  return indegree[v];
}

// Adding/removing edges

// adds an edge from v to w, validates the vertices, adds to the edge count
// adds w to v's adjacency list, and increases the in degree of w
void Digraph::add_edge(int v, int w) {
  validate_vertex(v);
  validate_vertex(w);
  _E++;
  _adj[v].push_back(w);
  indegree[w]++;
}

// removes an edge from v to w if it exists
// erases the index of w in v's adjacency list, decrements the number of edges, and subtracts w's in degree
void Digraph::remove_edge(int v, int w) {
  validate_vertex(v);
  validate_vertex(w);
  auto &v_list = this->_adj[v];
  auto it = std::find(v_list.begin(), v_list.end(), w);
  v_list.erase(it);
  this->_E--;
  indegree[w]--;
}

// creates a reverse graph of a given digraph (flips direction of all edges)
Digraph Digraph::reverse() const {
  Digraph r;
  r.V(BaseGraph::V());
  for (int v = 0; v < BaseGraph::V(); v++) {
    for (int w : this->adj(v)) {
      r.add_edge(w, v);
    }
  }

  return r;
}

// Cleanup
Digraph::~Digraph() noexcept{
  delete[] indegree;
}

/******************************************************************************
 *  Class: DepthFistSearch
 *  A class implementing the depth first search algorithm
 ******************************************************************************/
std::ostream& operator<<(std::ostream& out, Color &c){
  switch(c){
    case Color::White: return out << "W"; // White
    case Color::Grey:  return out << "G"; // Grey
    case Color::Black: return out << "B"; // Black
  }

  return out << "U"; // Unknown
}

// performs a depth first search on a BaseGraph g
// creates an array of vertex attributes of size v for the search
// vertex attributes: parent, color, time, and component
// for every vertex, checks that the color is white, then proceeds with DFS
DepthFirstSearch::DepthFirstSearch(BaseGraph &g): g(g), v_attributes(new VertexAttribute[g.V()]) {
  for(int v = 0; v < g.V(); v++){
    if(v_attributes[v].color == Color::White){
      dfs(v);
    }
  }
}

// same as above depth first search, but starts from a specific vertex
DepthFirstSearch::DepthFirstSearch(BaseGraph &g, int s): g(g), v_attributes(new VertexAttribute[g.V()]) {
  if (v_attributes[s].color == Color::White) {
    dfs(s);
  }
}

// same as above, except takes a list of source vertices to search from
DepthFirstSearch::DepthFirstSearch(BaseGraph& g, std::list<int> &sources): g(g), v_attributes(new VertexAttribute[g.V()]) {
  for(int s : sources) {
    if (v_attributes[s].color == Color::White) {
      dfs(s);
    }
  }
}

// performs a recursive dfs from vertex u
void DepthFirstSearch::dfs(int u) {
  // increments DFS time
  time++;
  // sets arrival time of u to the current time
  v_attributes[u].time[0] = time;
  // sets the color of u to grey (discovered but not done)
  v_attributes[u].color = Color::Grey;
  // if the graph g is directed, adds u to the preorder list
  if(g.is_directed()) pre.push_back(u);
  // for every vertex in u's adjacency list
  for(int v : g.adj(u)){
	  // if the vertex v is not discovered yet
    if(v_attributes[v].color == Color::White){
      // sets v's parent to u
      v_attributes[v].parent = u;
	  // performs a DFS from v
      dfs(v);
    }
  }
   // if graph g is directed, adds u to the postorder list (it is black at this point)
  if(g.is_directed()) post.push_back(u);
  // updates u's color to black
  v_attributes[u].color = Color::Black;
  // sets u's component to c_count
  v_attributes[u].component = c_count;
  // if vertex u has no parent, increments the component count
  // if there are other undiscovered vertices that haven't been searched from yet, ensures they will be in another component
  if(v_attributes[u].parent == -1){
    c_count++;
  }
  // increments class time
  time++;
 // sets u's end time to the current time
  v_attributes[u].time[1] = time;
}

// returns a stack of integers containing the path to vertex v from a component root vertex
std::stack<int> DepthFirstSearch::path_to(int v) {
  // creates an int stack to contain the path to vertex v
  std::stack<int> path;
  // sets the starting/current vertex to v
  int x = v;
  // while the current vertex has a parent
  while (v_attributes[x].parent != -1) {
	// add the current vertex to the stack
    path.push(x);
	// updates the current vertex
    x = v_attributes[x].parent;
  }
  // adds the final vertex (the parentless one) to the stack
  path.push(x);
  // returns the stack
  return path;
}

// returns the component containing vertex v
int DepthFirstSearch::component(int v){
  return v_attributes[v].component;
}

// returns the number of components
int DepthFirstSearch::components_count() {
  return c_count;
}

// returns whether vertex v was reached by the depth first search
bool DepthFirstSearch::reachable(int v) {
  return v_attributes[v].color == Color::Black;
}

// returns the preorder list compiled while performing the DFS
const std::list<int>& DepthFirstSearch::in_preorder() {
  return pre;
}

// returns the postorder list compiled while performing the DFS
const std::list<int>& DepthFirstSearch::in_postorder() {
  return post;
}

// returns the reverse of the postorder list
std::stack<int> DepthFirstSearch::in_reverse_postorder() {
  std::stack<int> reverse;
  for (int v : post) {
    reverse.push(v);
  }

  return reverse;
}

// prints the preorder list to the console
void DepthFirstSearch::show_in_preorder(std::ostream& out) {
  for (int v : pre) {
    out << v << " ";
  }

  out << std::endl;
}

// prints the postorder list to the console
void DepthFirstSearch::show_in_postorder(std::ostream& out) {
  for (int v : post) {
    out << v << " ";
  }
  
  out << std::endl;
}

// prints the reverse postorder list to the console
void DepthFirstSearch::show_in_reverse_postorder(std::ostream& out) {
  std::stack<int> reverse = in_reverse_postorder();
  while(!reverse.empty()){
    out << reverse.top() << " ";
    reverse.pop();
  }

  out << std::endl;
}

// string method for the search
std::string DepthFirstSearch::str()  const{
  std::ostringstream sout;
  for (int v = 0; v < g.V(); v++) {
    sout << v << ": (" << v_attributes[v].time[0]
         << '/' << v_attributes[v].time[1] << ")"
         << " - " << v_attributes[v].color
         << " - " << v_attributes[v].parent
         << " - " << v_attributes[v].component << '\n';
  }

  return sout.str();
}

DepthFirstSearch::~DepthFirstSearch() noexcept{
  delete[] v_attributes;
}