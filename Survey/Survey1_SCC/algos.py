
from collections import defaultdict
from collections import Counter
import queue


# Kosaraju's Algorithm Implementation
def kosaraju(graph):
    def dfs(graph, node, visited, stack=None):
        visited.add(node)
        for neighbor in graph[node]:
            if neighbor not in visited:
                dfs(graph, neighbor, visited, stack)
        if stack is not None:
            stack.append(node)

    # Step 1: Perform DFS and store nodes in finish time order
    stack = []
    visited = set()
    for node in graph:
        if node not in visited:
            dfs(graph, node, visited, stack)
    
    print(stack)
    # Step 2: Transpose the graph
    transpose = defaultdict(list)
    for node in graph:
        for neighbor in graph[node]:
            transpose[neighbor].append(node)
    
    # Step 3: Perform DFS on the transposed graph
    visited.clear()
    sccs = []
    
    while stack:
        node = stack.pop()
        if node not in visited:
            component = []
            dfs(transpose, node, visited, component)
            
            sccs.append(component)

    return sccs


def forward_backward_scc(graph):
    def bfs(start, graph):
        visited = set()
        queue = [start]
        visited.add(start)
        while queue:
            node = queue.pop(0)
            for neighbor in graph[node]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)
        return visited

    def reverse_graph(graph):
        reversed_graph = {node: [] for node in graph}
        for node in graph:
            for neighbor in graph[node]:
                reversed_graph[neighbor].append(node)
        return reversed_graph

    remaining_nodes = set(graph.keys())
    sccs = []
    reversed_graph = reverse_graph(graph)

    step = 1  # Step counter for debugging output
    while remaining_nodes:
        pivot = next(iter(remaining_nodes))
        print(f"\nStep {step}: Starting with pivot {pivot}")
        forward_reach = bfs(pivot, graph)
        print(f"Forward reachability from {pivot}: {forward_reach}")
        backward_reach = bfs(pivot, reversed_graph)
        print(f"Backward reachability to {pivot}: {backward_reach}")
        scc = forward_reach & backward_reach
        print(f"SCC identified with pivot {pivot}: {scc}")
        sccs.append(list(scc))
        remaining_nodes -= scc
        print(f"Remaining nodes after removing SCC: {remaining_nodes}")
        step += 1

    return sccs

def tarjan_scc(graph):
    def dfs(node):
        nonlocal time
        disc[node] = low[node] = time
        time += 1
        stack.append(node)
        stack_member[node] = True

        for neighbor in graph[node]:
            if disc[neighbor] == -1:  # If neighbor is not visited
                dfs(neighbor)
                low[node] = min(low[node], low[neighbor])
            elif stack_member[neighbor]:  # Back edge
                low[node] = min(low[node], disc[neighbor])

        # If node is a root of an SCC
        if low[node] == disc[node]:
            scc = []
            while stack:
                w = stack.pop()
                stack_member[w] = False
                scc.append(w)
                if w == node:
                    break
            sccs.append(scc)

    n = len(graph)
    disc = [-1] * n
    low = [-1] * n
    stack = []
    stack_member = [False] * n
    sccs = []
    time = 0

    for i in range(n):
        if disc[i] == -1:
            dfs(i)

    return sccs

def chain_algorithm(graph, K=None):
    """
    Chain Algorithm for finding Strongly Connected Components (SCCs),
    rewritten using clear definitions of `Post` and `Pre`, with step-by-step output.
    """

    def post(v, graph):
        """Return the direct successors of vertex `v`."""
        return set(graph[v])

    def pre(v, reversed_graph):
        """Return the direct predecessors of vertex `v`."""
        return set(reversed_graph[v])

    def reverse_graph(graph):
        """Reverse the edges of the graph."""
        reversed_graph = defaultdict(list)
        for node in graph:
            for neighbor in graph[node]:
                reversed_graph[neighbor].append(node)
        return reversed_graph

    def chain_recursive(graph, vertices, K, step=1):
        """Recursive function for finding SCCs."""
        if not vertices:
            #print(f"Step {step}: No vertices left to process. Returning.")
            return []

        #print(f"Step {step}: Current graph vertices: {vertices}")

        # Choose a pivot
        pivot = next(iter(K or vertices))
        #print(f"Step {step}: Chosen pivot: {pivot}")

        # Initialize forward and backward reachability sets
        forward_set = set()
        backward_set = set()
        layer = {pivot}
        reversed_graph = reverse_graph(graph)

        # Compute forward reachability (Fwd)
        #print(f"Step {step}: Computing forward reachability from pivot {pivot}")
        while layer:
            #print(f"Step {step}: Current forward layer: {layer}")
            forward_set.update(layer)
            new_layer = set()
            # for v in layer:
            #     new_layer.update(post(v, graph))
            
            # layer = new_layer - forward_set
            
            for v in layer:
                for neighbor in post(v, graph):
                    if neighbor not in forward_set:  # Avoid reprocessing
                        new_layer.add(neighbor)
            layer = new_layer
        
            print(f"Step {step}: Current forward layer: {layer}")
            
        # bfs
        forward_set_1 = set()
        q = queue.Queue()
        q.put(pivot)
        while not q.empty():
            node = q.get()
            forward_set_1.add(node)
            for neighbor in graph[node]:
                if neighbor not in forward_set_1:
                    q.put(neighbor)
                    
        # dfs
        forward_set_2 = set()
        stack = [pivot]
        while stack:
            node = stack.pop()
            forward_set_2.add(node)
            for neighbor in graph[node]:
                if neighbor not in forward_set_2:
                    stack.append(neighbor)
        
        print(f"Step {step}: Forward reachability (Fwd): {forward_set}")
        print(f"Step {step}: Forward reachability (Fwd): {forward_set_1}")
        print(f"Step {step}: Forward reachability (Fwd): {forward_set_2}")
        #print(f"Step {step}: Forward reachability (Fwd): {forward_set}")

        # Compute backward reachability (Bwd)
        #print(f"Step {step}: Computing backward reachability to pivot {pivot}")
        layer = {pivot}
        while layer:
            #print(f"Step {step}: Current backward layer: {layer}")
            backward_set.update(layer)
            new_layer = set()
            for v in layer:
                 new_layer.update(pre(v, reversed_graph))
            layer = new_layer - backward_set

        print(f"Step {step}: Backward reachability (Bwd): {backward_set}")

        # The SCC for this pivot is the intersection of forward and backward sets
        scc = forward_set_2 & backward_set
        
        
        #print(f"Step {step}: SCC identified: {scc}")
        result = [list(scc)]

        # Remaining vertices for recursive calls
        For_minus_scc = forward_set_2 - scc

        V_minus_foward = vertices - forward_set
        # Recursive call on the forward set minus the SCC
        
        print(f"Step {step}: Recursive call on the forward set minus the SCC {For_minus_scc}")
        result += chain_recursive(
            {v: [n for n in graph[v] if n in For_minus_scc] for v in For_minus_scc},
            For_minus_scc,
            None,
            step + 1
        )
        print(f"Step {step}: Recursive call on the forward set minus the SCC {V_minus_foward}")
        # Recursive call on the rest of the graph excluding the forward set
        result += chain_recursive(
            {v: [n for n in graph[v] if n in V_minus_foward] for v in V_minus_foward},
            V_minus_foward,
            None,
            step + 2
        )
        print()

        return result

    # Initial call to the recursive function
    return chain_recursive(graph, set(graph.keys()), K or set())


# Example usage
graph = {
    
    3: [4],
    4: [5],
    5: [6, 7],
    6: [3],
    7: [],
    0: [1],
    1: [2, 6],
    2: [5, 0],
}
graph = {
    0: [],
    1: [2],
    2: [3],
    3: [1],  # SCC 1: 1 -> 2 -> 3 -> 1
    4: [5],
    5: [6],
    6: [4],  # SCC 2: 4 -> 5 -> 6 -> 4
    7: [8],
    8: [9],
    9: [7],  # SCC 3: 7 -> 8 -> 9 -> 7
    3: [4],  # Connect SCC 1 to SCC 2
    6: [7],  # Connect SCC 2 to SCC 3
    #9: [1]   # Connect SCC 3 to SCC 1
}
sccs = chain_algorithm(graph)
print("\nFinal SCCs:", sccs)

# tarjan_sccs = tarjan_scc(graph)
# print("\nTarjan's Algorithm SCCs:", tarjan_sccs)

# fwdbwd_sccs = forward_backward_scc(graph)
# print("\nForward-Backward Algorithm SCCs:", fwdbwd_sccs)

# kosaraju_sccs = kosaraju(graph)
# print("\nKosaraju's Algorithm SCCs:", kosaraju_sccs)


# %%



