
from collections import deque
def topological_sort_verbose(graph):
    # Initialize in-degree dictionary
    in_degree = {node: 0 for node in graph}
    for nodes in graph.values():
        for neighbor in nodes:
            in_degree[neighbor] += 1

    print("Initial in-degrees:", in_degree)

    # Queue for vertices with no incoming edges
    queue = deque([node for node in graph if in_degree[node] == 0])
    print("Initial queue:", list(queue))
    
    topo_order = []
    
    while queue:
        current = queue.popleft()
        topo_order.append(current)
        print(f"Processing node: {current}, Current topological order: {topo_order}")
        
        for neighbor in graph[current]:
            in_degree[neighbor] -= 1
            print(f"Decreased in-degree of {neighbor} to {in_degree[neighbor]}")
            if in_degree[neighbor] == 0:
                queue.append(neighbor)
                print(f"Added {neighbor} to queue")

        print("Queue after processing:", list(queue))

    # Check if all nodes are in topo_order
    if len(topo_order) == len(graph):
        print("Final topological order:", topo_order)
        return topo_order
    else:
        raise ValueError("Graph has a cycle, topological sort not possible.")

# Example graph
graph = {
    'A': ['B'],
    'B': ['C', 'E'],
    'C': [],
    'D': [],
    'E': ['D']
}


count = 0
def dfs(graph, node, visited, stack):
    visited.add(node)
    global count 
    start_time = count
    count += 1
    for neighbor in graph[node]:
        print(f"Visiting node {neighbor} from {node} with {graph[node]} is in{ neighbor in visited}")
        if neighbor not in visited:
            dfs(graph, neighbor, visited, stack)
    stack.append([node, start_time])


# Excercise 1
graph = {
    1: [2],
    2:[3],
    3: [1, 4],
    4: [5],
    5: [6],
    6: [7],
    7: [5],
    8: [9],
    9: [10],
    10: [],
    11: [10],
    12: [11],
    12: [13],
    13: [14],
    14: []
    
}

# Excercise 2
graph ={
    'A': ['B', 'C'],
    'B': ['A'],
    'C': ['A']
}


graph = {
    'u':[],
    'v': ['u'],
    'w': ['v'],
}

visited = set()
stacks = []
count = 0
for node in graph:
    stack = []
    if node not in visited:
        dfs(graph, node, visited, stack)
        print(stack)
        stacks.append(stack)
print(stacks[::-1])



