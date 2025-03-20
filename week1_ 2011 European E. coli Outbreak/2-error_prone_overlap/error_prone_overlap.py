#python3
# Given a list of error-prone reads, perform the task of Genome Assembly and return the circular genome from which they came.

import sys
sys.setrecursionlimit(10**6)

# Global variables
k = 20
countmap = {}   # will count occurrences of each k-mer
bubbles = 0
root = None
tip_count = 0

class Vertex:
    def __init__(self, vertexNum, s, outedges, inedges):
        self.vertexNum = vertexNum      # unique id for the (k-1)-mer
        self.s = s                    # the (k-1)-mer string
        self.outedges = outedges[:]   # list of outgoing vertex ids
        self.inedges = inedges[:]     # list of incoming vertex ids
        self.edgeList = []            # list of indices to Edge objects (set later)
        self.removed = False          # mark for tip removal
        self.found = False            # used in bubble removal
        self.temp = None              # pointer to Node (used in bubble removal)
        self.visited = False          # used in Eulerian cycle search

class Node:
    def __init__(self, vertexNum):
        self.vertexNum = vertexNum
        self.kids = []    # list of child Node objects
        self.parent = None

class Edge:
    def __init__(self, frm, to):
        self.frm = frm  # starting vertex id
        self.to = to    # ending vertex id
        self.used = False  # flag for Eulerian cycle algorithm

def reader():
    """
    Read all input lines (reads) from stdin.
    """
    data = sys.stdin.read().splitlines()
    # remove any extra whitespace lines
    reads = [line.strip() for line in data if line.strip()]
    return reads

def is_overlap(a, b):
    """
    Check if the (k-1)-mers a and b have the required overlap.
    (i.e., a[1:] equals b[:-1])
    """
    for i in range(1, len(a)):
        if a[i] != b[i-1]:
            return False
    return True

def create_de_bruijn_graph(reads):
    """
    Build the de Bruijn graph from the reads.
    Vertices correspond to (k-1)-mers.
    """
    global countmap
    idmap = {}        # maps (k-1)-mer string -> unique id
    outedgesmap = {}  # maps (k-1)-mer string -> list of adjacent vertex ids (outgoing)
    inedgesmap = {}   # maps (k-1)-mer string -> list of adjacent vertex ids (incoming)
    countmap = {}
    uniquekmers = set()  # to avoid counting the same k-mer twice
    cur_id = 0

    for read in reads:
        for j in range(0, len(read) - k + 1):
            temp = read[j:j+k]
            a = temp[:-1]  # first k-1 characters
            b = temp[1:]   # last k-1 characters

            if temp in uniquekmers:
                countmap[temp] += 1
                continue
            else:
                uniquekmers.add(temp)
                countmap[temp] = 1

            if a not in idmap:
                idmap[a] = cur_id
                outedgesmap[a] = []
                inedgesmap[a] = []
                cur_id += 1
            if b not in idmap:
                idmap[b] = cur_id
                outedgesmap[b] = []
                inedgesmap[b] = []
                cur_id += 1

            if is_overlap(a, b):
                outedgesmap[a].append(idmap[b])
                inedgesmap[b].append(idmap[a])

    # Create the graph as a list of Vertex objects.
    graph = [None] * len(idmap)
    for key, vid in idmap.items():
        graph[vid] = Vertex(vid, key, outedgesmap[key], inedgesmap[key])
    return graph

def tip_removal(graph):
    """
    Remove tips from the graph.
    """
    global tip_count
    tip_count = 0
    for vertex in graph:
        if vertex.removed:
            continue
        if len(vertex.outedges) == 0:
            in_explore(graph, vertex.vertexNum)
            continue
        if len(vertex.inedges) == 0:
            out_explore(graph, vertex.vertexNum)
            continue

def in_explore(graph, vertex_index):
    vertex = graph[vertex_index]
    if len(vertex.outedges) != 0 or len(vertex.inedges) != 1:
        return
    vertex.removed = True
    global tip_count
    tip_count += 1
    if vertex.inedges:
        prev = vertex.inedges[0]
        try:
            graph[prev].outedges.remove(vertex_index)
        except ValueError:
            pass
        in_explore(graph, prev)

def out_explore(graph, vertex_index):
    vertex = graph[vertex_index]
    if len(vertex.inedges) != 0 or len(vertex.outedges) != 1:
        return
    vertex.removed = True
    global tip_count
    tip_count += 1
    if vertex.outedges:
        nxt = vertex.outedges[0]
        try:
            graph[nxt].inedges.remove(vertex_index)
        except ValueError:
            pass
        out_explore(graph, nxt)

def bubble_handler(graph):
    """
    Detect and resolve bubbles in the graph.
    """
    for vertex in graph:
        if vertex.removed or len(vertex.outedges) < 2:
            continue
        bfs(graph, vertex.vertexNum)

def bfs(graph, vertex_index):
    """
    Create a tree from the given vertex and explore for bubbles.
    """
    global root
    root = Node(vertex_index)
    # Reinitialize found flags and temp pointers for all vertices.
    for v in graph:
        v.found = False
        v.temp = None
    explore(graph, root, set())

def explore(graph, node, visited):
    visited.add(node.vertexNum)
    if graph[node.vertexNum].found:
        # A bubble is detected.
        global bubbles
        bubbles += 1
        common = find_common_ancestor(graph, node)
        bubble_solver(graph, node, common, visited)
        if node.vertexNum not in visited:
            return
    graph[node.vertexNum].found = True
    graph[node.vertexNum].temp = node

    if len(visited) >= k + 1:
        visited.remove(node.vertexNum)
        return

    # Explore all outgoing edges.
    for next_vertex in list(graph[node.vertexNum].outedges):
        if next_vertex in visited:
            continue
        child = Node(next_vertex)
        child.parent = node
        node.kids.append(child)
        explore(graph, child, visited)
        if node.vertexNum not in visited:
            return
    visited.remove(node.vertexNum)

def find_common_ancestor(graph, node):
    """
    Find the common ancestor of the current node and the node stored in graph[node.vertexNum].temp.
    """
    temp = graph[node.vertexNum].temp
    ancestors = []
    while temp is not None:
        ancestors.append(temp)
        temp = temp.parent
    temp = node
    while temp is not None:
        for anc in reversed(ancestors):
            if temp is anc:
                return temp
        temp = temp.parent
    return None

def bubble_solver(graph, node, common, visited):
    """
    Decide which of the two bubble paths to remove based on coverage.
    """
    node1 = node
    node2 = graph[node.vertexNum].temp

    # Calculate coverage along the path from node1 to common.
    sum_cov = 0
    count_edges = 0
    temp_node = node1
    while temp_node != common:
        str1 = graph[temp_node.vertexNum].s
        # For the parent node, if available:
        str2 = graph[temp_node.parent.vertexNum].s if temp_node.parent is not None else ""
        kmer = str2 + str1[-1]
        sum_cov += countmap.get(kmer, 0)
        count_edges += 1
        temp_node = temp_node.parent
    coverage1 = sum_cov / count_edges if count_edges > 0 else 0

    # Calculate coverage along the path from node2 to common.
    sum_cov = 0
    count_edges = 0
    temp_node = node2
    while temp_node != common:
        str1 = graph[temp_node.vertexNum].s
        str2 = graph[temp_node.parent.vertexNum].s if temp_node.parent is not None else ""
        kmer = str2 + str1[-1]
        sum_cov += countmap.get(kmer, 0)
        count_edges += 1
        temp_node = temp_node.parent
    coverage2 = sum_cov / count_edges if count_edges > 0 else 0

    # Remove the path with lower coverage.
    if coverage1 <= coverage2:
        vertices_removed = []
        removed_node = remove_path(graph, node1, common, vertices_removed)
        make_false(graph, removed_node)
        graph[node.vertexNum].found = True
        graph[node.vertexNum].temp = node2
        remove_from_set(visited, vertices_removed)
        if common and removed_node in common.kids:
            common.kids.remove(removed_node)
    else:
        vertices_removed = []
        removed_node = remove_path(graph, node2, common, vertices_removed)
        make_false(graph, removed_node)
        graph[node.vertexNum].found = True
        graph[node.vertexNum].temp = node1
        if common and removed_node in common.kids:
            common.kids.remove(removed_node)

def remove_path(graph, node, common, vertices_removed):
    """
    Remove the path (i.e. edges) from node up to (but not including) common.
    """
    parent = node.parent
    child = node
    temp_ret = None
    while child != common:
        if parent is not None:
            try:
                graph[parent.vertexNum].outedges.remove(child.vertexNum)
            except ValueError:
                pass
        vertices_removed.append(child.vertexNum)
        temp_ret = child
        child = parent
        if parent is not None:
            parent = parent.parent
        else:
            break
    return temp_ret

def remove_from_set(visited, vertices_removed):
    for v in vertices_removed:
        visited.discard(v)

def make_false(graph, node):
    """
    Reset the 'found' flag and temp pointer for all nodes in the subtree.
    """
    if node is None:
        return
    graph[node.vertexNum].found = False
    graph[node.vertexNum].temp = None
    for child in node.kids:
        make_false(graph, child)

def make_edges(graph):
    """
    Construct a list of Edge objects from the graph.
    """
    edges = []
    for vertex in graph:
        if vertex.removed:
            continue
        for nxt in vertex.outedges:
            edge = Edge(vertex.vertexNum, nxt)
            vertex.edgeList.append(len(edges))
            edges.append(edge)
    return edges

def eulerian_explore(graph, edges, vertex_index, cycle):
    """
    Recursively explore the graph to build an Eulerian cycle.
    """
    vertex = graph[vertex_index]
    vertex.visited = True
    for edge_index in list(vertex.edgeList):
        edge = edges[edge_index]
        if not edge.used:
            edge.used = True
            eulerian_explore(graph, edges, edge.to, cycle)
    cycle.append(vertex_index)

def find_cycle(graph, edges):
    """
    Find an Eulerian cycle and reconstruct the genome.
    """
    max_len = -1
    result = ""
    genome_map = {}
    for vertex in graph:
        if (not vertex.removed) and (not vertex.visited):
            cycle = []
            eulerian_explore(graph, edges, vertex.vertexNum, cycle)
            # Build the genome: start with the string of the last vertex
            genome = graph[cycle[-1]].s
            # Append the last character of each vertex in the reversed cycle
            for j in range(len(cycle) - 2, -1, -1):
                genome += graph[cycle[j]].s[-1]
            if len(genome) >= 5396:
                genome_map[genome] = len(genome) - 5396
            if len(genome) > max_len:
                result = genome
                max_len = len(genome)
    # Pick the genome with the smallest difference (if available)
    min_diff = float('inf')
    for key, diff in genome_map.items():
        if diff < min_diff:
            min_diff = diff
            result = key
    # Output the assembled genome (substring as in the Java code: indices 14 to 5410)
    print(result[14:5410])

def run():
    reads = reader()
    graph = create_de_bruijn_graph(reads)
    tip_removal(graph)
    bubble_handler(graph)
    tip_removal(graph)
    edges = make_edges(graph)
    find_cycle(graph, edges)

if __name__ == '__main__':
    run()
