# python3
# Genome Assembler for multiple datasets using de Bruijn graph approach

import sys

def parse_input():
  """Parses input. The first line is the number of reads (or read–pairs).
  Each subsequent line is either a single read or a read–pair (separated by '|').
  For paired reads, we simply add both reads to the list."""
  data = sys.stdin.read().splitlines()
  if not data:
    return []
  n = int(data[0])
  reads = []
  for line in data[1:]:
    line = line.strip()
    if not line:
      continue
    if '|' in line:
      parts = line.split('|')
      # parts[0] is READ1, parts[1] is READ2 (parts[2] is the expected gap, which we ignore)
      reads.append(parts[0])
      reads.append(parts[1])
    else:
      reads.append(line)
  return reads

def build_debruijn_graph(reads, k):
  """Builds a de Bruijn graph from all k-mers in the reads.
  Each k-mer produces an edge from its prefix (first k-1 letters) to its suffix (last k-1 letters).
  Returns:
   graph: dict mapping node->list of successor nodes
   in_degree: dict mapping node-># incoming edges
   out_degree: dict mapping node-># outgoing edges
  """
  graph = {}
  in_degree = {}
  out_degree = {}
  for read in reads:
    if len(read) < k:
      continue
    for i in range(len(read) - k + 1):
      kmer = read[i:i+k]
      prefix = kmer[:-1]
      suffix = kmer[1:]
      if prefix not in graph:
        graph[prefix] = []
      graph[prefix].append(suffix)
      out_degree[prefix] = out_degree.get(prefix, 0) + 1
      in_degree[suffix] = in_degree.get(suffix, 0) + 1
      # ensure both nodes are in the dictionaries and graph
      if suffix not in graph:
        graph[suffix] = []
      if prefix not in in_degree:
        in_degree[prefix] = in_degree.get(prefix, 0)
      if suffix not in out_degree:
        out_degree[suffix] = out_degree.get(suffix, 0)
  return graph, in_degree, out_degree

def remove_tips(graph, in_degree, out_degree):
  """Performs a simple iterative tip removal.
  A tip is defined here as a node with no incoming or no outgoing edges (but not isolated).
  For each tip, we remove its outgoing edges and update the degrees of its neighbors.
  (This is a simplified error-correction step.)"""
  changed = True
  while changed:
    changed = False
    tips = []
    for node in list(graph.keys()):
      # if a node has either no incoming or no outgoing edges (and is not completely isolated)
      if (in_degree.get(node, 0) == 0 or out_degree.get(node, 0) == 0) and (in_degree.get(node, 0) + out_degree.get(node, 0) > 0):
        tips.append(node)
    for tip in tips:
      # remove all outgoing edges from tip and update the in-degrees of its successors
      for succ in graph[tip]:
        in_degree[succ] -= 1
      out_degree[tip] = 0
      graph[tip] = []
      changed = True

def get_contigs(graph, in_degree, out_degree):
  """Finds maximal nonbranching paths (contigs) in the de Bruijn graph.
  For each node that is not 1-in-1-out, each outgoing edge begins a contig that is extended
  until a node is reached that is not 1-in-1-out. Also, isolated cycles are handled.
  Returns a list of contig strings."""
  contigs = []
  def path_to_contig(path):
    contig = path[0]
    for node in path[1:]:
      contig += node[-1]
    return contig

  # Collect all nodes that appear as keys or as successors.
  nodes = set(graph.keys())
  for succs in graph.values():
    for v in succs:
      nodes.add(v)
   
  # Make a copy of the graph (we will remove edges as we traverse)
  graph_copy = {}
  for node, succs in graph.items():
    graph_copy[node] = list(succs)
   
  # For each node that is not 1-in-1, start contigs for each outgoing edge.
  for node in nodes:
    if out_degree.get(node, 0) > 0:
      if not (in_degree.get(node, 0) == 1 and out_degree.get(node, 0) == 1):
        while graph_copy[node]:
          next_node = graph_copy[node].pop(0)
          path = [node, next_node]
          current = next_node
          while in_degree.get(current, 0) == 1 and out_degree.get(current, 0) == 1:
            # current’s only successor
            if not graph_copy[current]:
              break
            next_current = graph_copy[current].pop(0)
            path.append(next_current)
            current = next_current
          contigs.append(path_to_contig(path))
  # Now, add isolated cycles if any remain.
  for node in nodes:
    if graph_copy.get(node, []):
      cycle = [node]
      current = node
      while True:
        if not graph_copy[current]:
          break
        next_node = graph_copy[current].pop(0)
        cycle.append(next_node)
        current = next_node
        if current == node:
          break
      if len(cycle) > 1:
        contigs.append(path_to_contig(cycle))
  return contigs

def main():
  reads = parse_input()
  if not reads:
    return
  # For simplicity, choose k equal to the length of the first read.
  # (In many datasets the read lengths are uniform; for paired reads, we assume both have similar length.)
  k = len(reads[0])
  for r in reads:
    if len(r) < k:
      k = len(r)
  graph, in_degree, out_degree = build_debruijn_graph(reads, k)
  # Apply tip removal as an error-correction step.
  remove_tips(graph, in_degree, out_degree)
  contigs = get_contigs(graph, in_degree, out_degree)
  # Output contigs in FASTA format.
  for i, contig in enumerate(contigs):
    print(">CONTIG{}".format(i+1))
    print(contig)

if __name__ == "__main__":
  main()