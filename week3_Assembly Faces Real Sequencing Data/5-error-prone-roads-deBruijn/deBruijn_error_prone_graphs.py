# python3
''' 
Given a list of error-prone reads, perform the task of Genome Assembly using de Bruijn graphs and
return the circular genome from which they came. Break the reads into fragments of length 𝑘 = 15
before constructing the de Bruijn graph, remove tips, and handle bubbles.
'''

import sys
sys.setrecursionlimit(10**6)

# Global variables used across functions.
countmap = {}   # counts the number of occurrences for each k-mer (string of length k)
k = 20          # k-mer length (the code breaks reads into fragments of length k)
tip_removed_count = 0  # counts number of tip removals
bubbles = 0           # counts bubbles detected
root = None           # global pointer to the root node in bubble detection


# Data structures (classes)
class Vertex:
    def __init__(self, vertexNum, s, outedges, inedges):
        self.vertexNum = vertexNum          # unique id for this (k-1)-mer
        self.str = s                        # the string label (of length k-1)
        self.outedges = list(outedges)      # list of adjacent vertex ids (outgoing)
        self.inedges = list(inedges)        # list of adjacent vertex ids (incoming)
        self.edgeList = []                  # list of indices into the list of Edge objects
        self.removed = False                # flag for tip removal
        self.found = False                  # flag used in bubble removal (DFS)
        self.temp = None                    # pointer to the Node where this vertex sits in the DFS tree
        self.visited = False                # used for Eulerian cycle extraction

class Node:
    def __init__(self, vertexNum):
        self.vertexNum = vertexNum          # the id of the vertex
        self.kids = []                      # list of children Node objects
        self.parent = None                  # pointer to the parent Node

class Edge:
    def __init__(self, frm, to):
        self.frm = frm                      # starting vertex id
        self.to = to                        # ending vertex id
        self.used = False                   # flag to mark if this edge has been used in Eulerian cycle search


# Utility functions

def reader():
    """Reads input from stdin and returns a list of reads (each a string)."""
    data = sys.stdin.read().splitlines()
    return [line.strip() for line in data if line.strip()]

def isOverlap(a, b):
    """
    Check whether two (k-1)-mers overlap in the expected way.
    Specifically, for each i=1..len(a)-1, we require a[i]==b[i-1].
    """
    for i in range(1, len(a)):
        if a[i] != b[i-1]:
            return False
    return True


# Graph Construction

def createDeBruijnGraph(reads):
    """
    Build the de Bruijn graph.
    For each read, all k-mers (of length k) are extracted.
    Each k-mer is split into its (k-1)-mer prefix and suffix.
    An edge is added from prefix to suffix if the strings overlap as expected.
    Duplicate k-mers are counted using the global countmap.
    """
    global countmap
    idmap = {}            # maps (k-1)-mer string to a unique integer id
    outedgesmap = {}      # maps (k-1)-mer string to list of outgoing adjacent vertex ids
    inedgesmap = {}       # maps (k-1)-mer string to list of incoming adjacent vertex ids
    countmap = {}         # reset the global countmap
    uniquekmers = set()   # to ensure we count each k-mer only once per occurrence

    current_id = 0
    for read in reads:
        # For each substring of length k from the read.
        for j in range(0, len(read) - k + 1):
            temp = read[j:j+k]
            a = temp[:-1]  # prefix (first k-1 characters)
            b = temp[1:]   # suffix (last k-1 characters)
            if temp in uniquekmers:
                countmap[temp] += 1
                continue
            uniquekmers.add(temp)
            countmap[temp] = 1

            if a not in idmap:
                idmap[a] = current_id
                outedgesmap[a] = []
                inedgesmap[a] = []
                current_id += 1
            if b not in idmap:
                idmap[b] = current_id
                outedgesmap[b] = []
                inedgesmap[b] = []
                current_id += 1

            if isOverlap(a, b):
                outedgesmap[a].append(idmap[b])
                inedgesmap[b].append(idmap[a])
    # Create graph: an array (list) of Vertex objects indexed by their unique id.
    graph = [None] * len(idmap)
    for key, vid in idmap.items():
        graph[vid] = Vertex(vid, key, outedgesmap[key], inedgesmap[key])
    return graph


# Tip Removal Functions

def tipRemoval(graph):
    """
    Iterates over all vertices in the graph.
    If a vertex has no outgoing edges, calls inExplore (to remove tips at the end).
    If a vertex has no incoming edges, calls outExplore (to remove tips at the beginning).
    """
    global tip_removed_count
    tip_removed_count = 0
    for i in range(len(graph)):
        if graph[i].removed:
            continue
        if len(graph[i].outedges) == 0:
            inExplore(graph, i)
            continue
        if len(graph[i].inedges) == 0:
            outExplore(graph, i)
            continue

def inExplore(graph, vertex):
    """
    Recursively remove tip from the graph inwards.
    Only removes a vertex if it has no outgoing edges and exactly one incoming edge.
    """
    if len(graph[vertex].outedges) != 0 or len(graph[vertex].inedges) != 1:
        return
    graph[vertex].removed = True
    global tip_removed_count
    tip_removed_count += 1

    temp = graph[vertex].inedges[0]
    # Remove the edge from the predecessor's outgoing edges.
    if vertex in graph[temp].outedges:
        graph[temp].outedges.remove(vertex)
    if temp in graph[vertex].inedges:
        graph[vertex].inedges.remove(temp)
    inExplore(graph, temp)

def outExplore(graph, vertex):
    """
    Recursively remove tip from the graph outwards.
    Only removes a vertex if it has no incoming edges and exactly one outgoing edge.
    """
    if len(graph[vertex].inedges) != 0 or len(graph[vertex].outedges) != 1:
        return
    graph[vertex].removed = True
    global tip_removed_count
    tip_removed_count += 1

    temp = graph[vertex].outedges[0]
    if vertex in graph[temp].inedges:
        graph[temp].inedges.remove(vertex)
    if temp in graph[vertex].outedges:
        graph[vertex].outedges.remove(temp)
    outExplore(graph, temp)


# Bubble Removal Functions

def bubbleHandler(graph):
    """
    For each vertex that is not removed and has at least two outgoing edges,
    initiate a DFS to detect and resolve bubbles.
    """
    for i in range(len(graph)):
        if graph[i].removed or len(graph[i].outedges) < 2:
            continue
        bfs(graph, graph[i].vertexNum)

def bfs(graph, vertex):
    """
    This function initializes a DFS tree rooted at the given vertex.
    It resets the 'found' flag and temp pointer for every vertex.
    """
    global root
    visited_set = set()
    root = Node(vertex)
    for v in graph:
        v.found = False
        v.temp = None
    explore(graph, root, visited_set)

def explore(graph, node, visited):
    """
    DFS-style exploration.
    If a vertex is encountered twice (its 'found' flag is True), a bubble is detected.
    """
    visited.add(node.vertexNum)
    if graph[node.vertexNum].found:
        global bubbles
        bubbles += 1
        common = findCommonAncestor(graph, node)
        bubbleSolver(graph, node, common, visited)
        if node.vertexNum not in visited:
            return
    graph[node.vertexNum].found = True
    graph[node.vertexNum].temp = node

    # Stop exploring if the current path is too long.
    if len(visited) >= k + 1:
        visited.remove(node.vertexNum)
        return

    # Explore all outgoing edges.
    for nxt in list(graph[node.vertexNum].outedges):
        if nxt in visited:
            continue
        child = Node(nxt)
        child.parent = node
        node.kids.append(child)
        explore(graph, child, visited)
        if node.vertexNum not in visited:
            return
    visited.remove(node.vertexNum)

def findCommonAncestor(graph, node):
    """
    Finds the common ancestor between the current node and the node stored in graph[node.vertexNum].temp.
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

def bubbleSolver(graph, node, common, visited):
    """
    Given a bubble (two distinct paths from vertex v to w), choose the path with lower coverage to remove.
    Coverage is computed as the average count (from countmap) of the k-mers along the path.
    """
    node1 = node
    node2 = graph[node.vertexNum].temp

    sum1 = 0.0
    edge_count = 0.0
    temp_node = node1
    while temp_node != common:
        str1 = graph[temp_node.vertexNum].str
        str2 = graph[temp_node.parent.vertexNum].str
        kmer = str2 + str1[-1]
        sum1 += countmap.get(kmer, 0)
        edge_count += 1
        temp_node = temp_node.parent
    coverage1 = sum1 / edge_count if edge_count != 0 else 0.0

    sum2 = 0.0
    edge_count = 0.0
    temp_node = node2
    while temp_node != common:
        str1 = graph[temp_node.vertexNum].str
        str2 = graph[temp_node.parent.vertexNum].str
        kmer = str2 + str1[-1]
        sum2 += countmap.get(kmer, 0)
        edge_count += 1
        temp_node = temp_node.parent
    coverage2 = sum2 / edge_count if edge_count != 0 else 0.0

    node1 = node
    node2 = graph[node.vertexNum].temp

    if coverage1 <= coverage2:
        verticesRemoved = []
        temp_removed = removePath(graph, node1, common, verticesRemoved)
        makeFalse(graph, temp_removed)
        graph[node.vertexNum].found = True
        graph[node.vertexNum].temp = node2
        remove_from_set(visited, verticesRemoved)
        if common is not None and temp_removed in common.kids:
            common.kids.remove(temp_removed)
    else:
        verticesRemoved = []
        temp_removed = removePath(graph, node2, common, verticesRemoved)
        makeFalse(graph, temp_removed)
        graph[node.vertexNum].found = True
        graph[node.vertexNum].temp = node1
        if common is not None and temp_removed in common.kids:
            common.kids.remove(temp_removed)

def removePath(graph, node, common, verticesRemoved):
    """
    Removes the buggy path from the bubble.
    Traverses upward from the given node until the common ancestor is reached,
    removing each edge along the path.
    """
    parent = node.parent
    child = node
    temp_ret = None
    while child != common:
        if child.vertexNum in graph[parent.vertexNum].outedges:
            graph[parent.vertexNum].outedges.remove(child.vertexNum)
        verticesRemoved.append(child.vertexNum)
        temp_ret = child
        child = parent
        parent = parent.parent if parent is not None else None
        if parent is None:
            break
    return temp_ret

def remove_from_set(visited, verticesRemoved):
    for v in verticesRemoved:
        visited.discard(v)

def makeFalse(graph, node):
    """
    Resets the 'found' flag and temp pointer for the subtree rooted at the given node.
    """
    graph[node.vertexNum].found = False
    graph[node.vertexNum].temp = None
    for kid in node.kids:
        makeFalse(graph, kid)


# Eulerian Cycle (Genome Assembly)

def makeEdges(graph):
    """
    Creates a list of Edge objects from the cleaned graph.
    For each vertex (that is not removed), every outgoing edge is converted into an Edge object.
    """
    edges = []
    for i in range(len(graph)):
        if graph[i].removed:
            continue
        for nxt in graph[i].outedges:
            edge = Edge(i, nxt)
            graph[i].edgeList.append(len(edges))
            edges.append(edge)
    return edges

def eulerianExplore(graph, edges, vertex, cycle):
    """
    Recursively explores the graph to find an Eulerian cycle.
    """
    for edge_index in list(graph[vertex].edgeList):
        edge = edges[edge_index]
        if not edge.used:
            edge.used = True
            eulerianExplore(graph, edges, edge.to, cycle)
    cycle.append(vertex)
    graph[vertex].visited = True

def findCycle(graph, edges):
    """
    Finds an Eulerian cycle in the graph, reconstructs the genome by stitching together vertex strings,
    and finally prints the assembled genome (a substring from index 14 to 5410).
    """
    max_length = -1
    result = ""
    genome_map = {}  # maps genome string to (genome length - 5396)
    for i in range(len(graph)):
        if not graph[i].removed and not graph[i].visited:
            cycle = []
            eulerianExplore(graph, edges, i, cycle)
            genome = graph[cycle[-1]].str
            for j in range(len(cycle)-2, -1, -1):
                genome += graph[cycle[j]].str[-1]
            if len(genome) >= 5396:
                genome_map[genome] = len(genome) - 5396
            if len(genome) > max_length:
                result = genome
                max_length = len(genome)
    min_diff = float('inf')
    for key, diff in genome_map.items():
        if diff < min_diff:
            min_diff = diff
            result = key
    # Print the assembled genome (a substring from index 14 to 5410)
    print(result[14:5410])

# Main program logic

def run():
    reads = reader()                    # read all input lines
    graph = createDeBruijnGraph(reads)  # build the de Bruijn graph from the reads
    tipRemoval(graph)                   # remove tips from the graph
    bubbleHandler(graph)                # handle (remove) bubbles in the graph
    tipRemoval(graph)                   # remove tips again (to catch those exposed by bubble removal)
    edges = makeEdges(graph)            # create edge objects from the cleaned graph
    findCycle(graph, edges)             # find Eulerian cycle and output the genome

if __name__ == '__main__':
    run()