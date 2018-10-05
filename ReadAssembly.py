class AssemblyGraph:
    def __init__(self, vertices=None):
        """Constructs a directed graph with num_vertices vertices and zero edges"""
        self.vertices = dict() if not vertices else self.add_vertices(vertices)
        self.edges = dict()
        self.connectivity = 0  # keep track of this for efficiency

    def add_vertex(self, vertex):
        """obligatory add_vertex method"""
        self.vertices[vertex] = {"in": [], "out": []}

    def add_vertices(self, vertices):
        """adds all vertices in a list to our set"""
        self.vertices = dict((v, {"in": [], "out": []}) for v in vertices)
        return self.vertices  # return the list we create so that there's a return

    def get_edge(self, i, j):
        """Returns edge weight if the graph contains the directed edge (i, j), null otherwise."""
        return self.edges.get("{} {}".format(i, j))

    @staticmethod
    def overlap_length(left, right):
        """Returns the length of the longest suffix of left that is a prefix of right"""
        ret = 0
        for i in range(1, min(len(left), len(right)) + 1):
            if left.endswith(right[:i]):
                ret = i
        return ret

    def add_edge(self, i, j):
        """Adds the directed edge (i, j) to the graph."""
        self.edges["{} {}".format(i, j)] = self.overlap_length(i, j)
        self.vertices[i]["out"].append(j)
        self.vertices[j]["in"].append(i)
        self.connectivity += 2

    def out_edges(self, i):
        """Returns a list of directed edges outgoing from vertex i."""
        return [(i, j) for j in self.vertices[i]["out"]]

    def in_edges(self, j):
        """Returns a list of directed edges incoming to vertex j."""
        return [(i, j) for i in self.vertices[j]["in"]]

    def outdegree(self, i):
        """Returns the outdegree of vertex i."""
        return len(self.vertices[i]["out"])

    def indegree(self, i):
        """Returns the indegree of vertex i."""
        return len(self.vertices[i]["in"])

    def degree(self, i):
        """Returns the degree of vertex i."""
        return self.indegree(i) + self.outdegree(i)

    def add_edges(self, edges):
        """Adds all edges from a list to the graph."""
        for i, j in edges:
            self.add_edge(i, j)

    def num_vertices(self):
        """Returns the number of vertices in the graph."""
        return len(self.vertices)

    def num_edges(self):
        """Returns the number of edges in the graph."""
        return len(self.edges)

    def get_edges(self):
        """Returns an iterator over the edges of the graph."""
        return [edge.split(" ") for edge in self.edges.keys()]

    def generate_possible_edges(self):
        """optimizes sorting after generating all possible edges"""
        from operator import itemgetter
        edges = [(i, j, self.overlap_length(i, j)) for j in self.vertices for i in self.vertices if i != j]
        lexi_sort = sorted(edges, key=itemgetter(0, 1), reverse=True)
        weight_sort = sorted(lexi_sort, key=itemgetter(2))
        return weight_sort

    def is_disconnected(self):
        """uses math to determine if the directed graph is connected or not
        """
        degrees_per_edge = 2
        num_edges_in_connected_graph = (self.num_vertices() - 1)
        full_connectivity = num_edges_in_connected_graph * degrees_per_edge
        return self.connectivity != full_connectivity

    def does_not_cause_cycle(self, edge):
        """starts at an edge's destination and then looks through all forward-connected nodes
        to make sure that the edge source isn't present
        """
        src, dest = edge
        to_check = list(self.vertices[dest]["out"])  # look up out destinations (MAKE A COPY)
        while len(to_check) != 0:
            curr = to_check.pop()
            if curr == src:  # we found a cycle!
                return False
            to_check.extend(self.vertices[curr]["out"])  # add these destinations
        return True

    def merge_ordered_reads(self, reads):
        """Returns the shortest superstring resulting from merging
        the elements of reads, an ordered list of strings"""
        # create list of overlaps
        overlaps = [self.overlap_length(reads[i], reads[i + 1]) for i in range(len(reads) - 1)]
        # merge
        superstring = ""
        for i, read in enumerate(reads):
            if i == 0:  # base case
                superstring = read
                continue

            superstring += read[overlaps[i - 1]:]  # reference the overlap list
        return superstring

    def superstring_from_edges(self):
        """ figures out which of the reads are at the tip of the superstring, then determines
         which tip is which. after that, it just sequentially reassembles and then merges"""
        from collections import Counter
        from functools import reduce

        # figure out which end is which in our edges
        duplicate_list = reduce(lambda x, y: x + y, self.get_edges())  # unpack all the edge tuples
        frequencies = dict(Counter(duplicate_list)).items()  # count up frequencies
        superstring_tips = filter(lambda x: x[1] == 1, frequencies)  # grab the two reads that only occur once
        source, terminal = (t[0] for t in superstring_tips)  # arbitrarily assign them

        # check if we got it backwards and swap if so
        sources, terminals = zip(*self.get_edges())
        if source not in sources:
            source, terminal = terminal, source

        # assemble!
        read_mapper = dict((k, v) for k, v in self.get_edges())

        curr = source
        ordered_reads = [curr]
        while curr != terminal: # stop once we've attached the terminal end 
            ordered_reads.append(read_mapper[curr])
            curr = read_mapper[curr]
        return self.merge_ordered_reads(ordered_reads)

    
from collections import deque as queue

def greedy_assemble(reads):
    """Returns a string that is a superstring of the input reads, which are given as a list of strings.
    The superstring computed is determined by the greedy algorithm as described in HW1, with specific tie-breaking
    criteria.
    """
    g = AssemblyGraph(reads)
    q = queue(g.generate_possible_edges())
    while g.is_disconnected():
        src, dest = edge = q.pop()[0:2]  # pull out the two vertices, ignore weight
        if g.outdegree(src) == g.indegree(dest) == 0 and g.does_not_cause_cycle(edge):
            g.add_edge(*edge)

    return g.superstring_from_edges()
