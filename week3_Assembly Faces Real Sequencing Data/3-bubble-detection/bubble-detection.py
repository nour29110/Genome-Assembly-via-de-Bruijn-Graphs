# python3

''' 
Given a list of error-prone reads and two integers, 𝑘 and 𝑡, construct a de Bruijn graph from the
𝑘-mers created from the reads and perform the task of bubble detection on this de Bruijn graph with
a path length threshold of 𝑡.
'''

import sys

class DeBrujin:
    def __init__(self, k, t):
        self.k = k                   # k-mer size
        self.t = t                   # bubble length threshold
        self.v = 0                   # counter for vertex ids
        self.edges = []              # list of unique k-mers (edges)
        self.edge_set = set()        # set to quickly check if a k-mer has been seen
        self.vertices = []           # list of vertex strings (k-1-mers)
        self.vertex_map = {}         # maps vertex string to its id (index in vertices)
        self.adjList = []            # adjacency list: list of lists of ints

    def add_read(self, r):
        # For each read, extract every k-mer and add if not already seen.
        for i in range(len(r) - self.k + 1):
            kmer = r[i:i + self.k]
            if kmer not in self.edge_set:
                self.edge_set.add(kmer)
                self.edges.append(kmer)

    def build_de_brujin(self):
        # For every edge (k-mer), break it into prefix and suffix (each of length k-1)
        for edge in self.edges:
            pre = edge[:-1]
            suf = edge[1:]
            if pre not in self.vertex_map:
                self.vertices.append(pre)
                self.vertex_map[pre] = self.v
                self.v += 1
            if suf not in self.vertex_map:
                self.vertices.append(suf)
                self.vertex_map[suf] = self.v
                self.v += 1
        # Initialize the adjacency list for all vertices.
        self.adjList = [[] for _ in range(self.v)]
        # For each edge, add a directed edge from its prefix vertex to its suffix vertex.
        for edge in self.edges:
            pre = edge[:-1]
            suf = edge[1:]
            from_vertex = self.vertex_map[pre]
            to_vertex = self.vertex_map[suf]
            self.adjList[from_vertex].append(to_vertex)

    def count_bubbles(self):
        cnt = 0
        # For every vertex with at least 2 outgoing edges,
        # count bubbles starting from that source.
        for i in range(len(self.adjList)):
            if len(self.adjList[i]) >= 2:
                cnt += self.count_bubbles_from_source(i, self.t)
        return cnt

    def count_bubbles_from_source(self, s, t):
        count = 0
        nPath = len(self.adjList[s])
        # For every pair of distinct outgoing edges from s:
        for i in range(nPath - 1):
            for j in range(i + 1, nPath):
                left_paths = []
                left_sets = []
                visited = set()
                path = [self.adjList[s][i]]
                visited.add(self.adjList[s][i])
                self.get_non_overlapping_paths(path, visited, left_paths, left_sets, t)
                
                right_paths = []
                right_sets = []
                visited2 = set()
                path2 = [self.adjList[s][j]]
                visited2.add(self.adjList[s][j])
                self.get_non_overlapping_paths(path2, visited2, right_paths, right_sets, t)
                
                count += self.count_bubbles_from_left_right_paths(left_sets, right_sets, left_paths, right_paths)
        return count

    def get_non_overlapping_paths(self, path, visited, all_paths, all_sets, l):
        # If the current path has reached length l, record it.
        if len(path) == l:
            all_paths.append(path.copy())
            all_sets.append(visited.copy())
            return
        s = path[-1]
        # If the current vertex has no outgoing edge, record the path.
        if len(self.adjList[s]) == 0:
            all_paths.append(path.copy())
            all_sets.append(visited.copy())
            return
        # Extend the path for each unvisited neighbor.
        for v in self.adjList[s]:
            if v not in visited:
                visited.add(v)
                path.append(v)
                self.get_non_overlapping_paths(path, visited, all_paths, all_sets, l)
                path.pop()
                visited.remove(v)
            else:
                continue

    def get_path_string(self, left_path, right_path, merge):
        ret = ""
        for p in left_path:
            if p == merge:
                ret += str(p)
                break
            else:
                ret += str(p) + ","
        ret += ";"
        for p in right_path:
            if p == merge:
                ret += str(p)
                break
            else:
                ret += str(p) + ","
        return ret

    def count_bubbles_from_left_right_paths(self, left_sets, right_sets, left_paths, right_paths):
        merge_paths = set()
        # For each left path and each right path, check for a common vertex (the merge point).
        for i, left_set in enumerate(left_sets):
            for right_path in right_paths:
                merged = False
                for v in right_path:
                    if v in left_set:
                        merged = True
                        path_string = self.get_path_string(left_paths[i], right_path, v)
                        merge_paths.add(path_string)
                        break  # Found a merge point for this pair; go to the next right path.
                if merged:
                    continue
        return len(merge_paths)

def main():
    data = sys.stdin.read().splitlines()
    if not data:
        return
    # The first line contains k and t separated by a space.
    parts = data[0].split()
    k = int(parts[0])
    t = int(parts[1])
    graph = DeBrujin(k, t)
    # Process each subsequent line as a read.
    for line in data[1:]:
        read = line.strip()
        if read:
            graph.add_read(read)
    graph.build_de_brujin()
    result = graph.count_bubbles()
    print(result)

if __name__ == '__main__':
    main()