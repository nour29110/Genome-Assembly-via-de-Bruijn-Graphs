# python3
# Given a list of error-prone reads, construct a de Bruijn graph from the 15-mers created from the reads and perform the task of tip removal on this de Bruijn graph.

import sys
from collections import defaultdict

sys.setrecursionlimit(10**7)

class Vertex:
    def __init__(self, label):
        self.label = label         # (k-1)-mer label
        self.inedges = []          # list of incoming vertex indices
        self.outedges = []         # list of outgoing vertex indices
        self.removed = False       # flag to indicate if this vertex is removed

def build_debruijn_graph(reads, k):
    """
    Build a de Bruijn graph from all k-mers extracted from the reads.
    Each k-mer is split into prefix (k-1) and suffix (k-1),
    and an edge is added from prefix -> suffix.
    Returns:
      graph: list of Vertex objects
      mapping: dict { (k-1)-mer_str -> vertex_index }
    """
    mapping = {}
    graph = []
    def get_vertex(label):
        """Return the index of the vertex with given label, creating if needed."""
        if label in mapping:
            return mapping[label]
        else:
            idx = len(graph)
            graph.append(Vertex(label))
            mapping[label] = idx
            return idx

    for read in reads:
        read = read.strip()
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            prefix = kmer[:-1]
            suffix = kmer[1:]
            v1 = get_vertex(prefix)
            v2 = get_vertex(suffix)
            graph[v1].outedges.append(v2)
            graph[v2].inedges.append(v1)

    return graph, mapping

def in_explore(graph, idx):
    """
    Recursively remove a 'right tip' (vertex with no outgoing edges, exactly 1 incoming edge).
    Returns the number of edges removed in this chain.
    """
    v = graph[idx]
    if v.removed or len(v.outedges) != 0 or len(v.inedges) != 1:
        return 0
    # Mark the vertex as removed
    v.removed = True
    # There's exactly one incoming edge
    pred_id = v.inedges[0]
    pred = graph[pred_id]
    # Remove the edge from pred -> v
    removed_edges = 0
    if idx in pred.outedges:
        pred.outedges.remove(idx)
        removed_edges += 1
    # Remove the back-pointer
    v.inedges.remove(pred_id)
    # Recurse on the predecessor
    return removed_edges + in_explore(graph, pred_id)

def out_explore(graph, idx):
    """
    Recursively remove a 'left tip' (vertex with no incoming edges, exactly 1 outgoing edge).
    Returns the number of edges removed in this chain.
    """
    v = graph[idx]
    if v.removed or len(v.inedges) != 0 or len(v.outedges) != 1:
        return 0
    # Mark the vertex as removed
    v.removed = True
    # There's exactly one outgoing edge
    succ_id = v.outedges[0]
    succ = graph[succ_id]
    # Remove the edge from v -> succ
    removed_edges = 0
    if idx in succ.inedges:
        succ.inedges.remove(idx)
        removed_edges += 1
    v.outedges.remove(succ_id)
    # Recurse on the successor
    return removed_edges + out_explore(graph, succ_id)

def tip_removal(graph):
    """
    Iteratively remove all tips until no more tips can be removed in a pass.
    Returns the total number of edges removed.
    """
    total_removed = 0
    while True:
        changed = False
        # Iterate over a snapshot of vertices so we don't skip newly exposed tips
        for i, v in enumerate(graph):
            if v.removed:
                continue
            # If no outedges => possible right tip
            if len(v.outedges) == 0:
                edges_removed = in_explore(graph, i)
                if edges_removed > 0:
                    total_removed += edges_removed
                    changed = True
            # Else if no inedges => possible left tip
            elif len(v.inedges) == 0:
                edges_removed = out_explore(graph, i)
                if edges_removed > 0:
                    total_removed += edges_removed
                    changed = True
        if not changed:
            break
    return total_removed

def main():
    # Read the 400 lines (each is a 100-nt read).
    reads = sys.stdin.read().split()
    # k = 15 (as stated in the problem)
    k = 15
    # Build the de Bruijn graph
    graph, _ = build_debruijn_graph(reads, k)
    # Perform iterative tip removal
    removed_count = tip_removal(graph)
    # Print the final number of edges removed
    print(removed_count)

if __name__ == '__main__':
    main()
