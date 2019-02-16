import logging
from pygssw.gssw import gssw_node_create, gssw_create_nt_table, gssw_create_score_matrix, \
    gssw_graph_create, gssw_create_score_matrix, gssw_nodes_add_edge, gssw_graph_add_node, \
    gssw_graph_fill, gssw_graph_print_score_matrices, gssw_graph_trace_back, test_wrapper


class Aligner:
    def __init__(self, ob_graph, sequence_graph, start_node, sequence):
        self.ob_graph = ob_graph
        self.start_node = start_node
        self._sequence_graph = sequence_graph
        self.sequence = sequence
        self._n_bp_to_traverse = len(sequence) + 10

        self.is_reverse = False
        self.adj_list = ob_graph.adj_list
        if self.start_node < 0:
            self.is_reverse = True
            self.adj_lis = ob_graph.reverse_adj_list


        self._nodes = set()
        self._edges = set()

        self._traversed_nodes = {}  # node id => bp to that node

        self._find_all_local_nodes_and_edges()

    def _find_all_local_nodes_and_edges(self):
        self._traverse_from_node(self.start_node, None, n_bp_traversed=0)

    def _traverse_from_node(self, node, prev_node, n_bp_traversed):
        #print("On node %d. Prev: %s" % (node, prev_node))
        node = int(node)
        if prev_node is not None:
            #print("Adding edge %d-%d" % (prev_node, node))
            if node < prev_node:
                logging.error("Graph is not sorted. Only sorted graphs currently supported")

            self._edges.add((prev_node, node))

        if node in self._traversed_nodes:  # and self._traversed_nodes[node] <= n_bp_traversed:
            #print("   Stopping. already processed")
            # Already processed this node with a possibly shorter path
            return

        self._traversed_nodes[node] = n_bp_traversed

        self._nodes.add(node)

        n_bp_traversed += self.ob_graph.blocks[node].length()

        if n_bp_traversed > self._n_bp_to_traverse:
            #print("    Stopping because length")
            return

        edges_out = self.adj_list[node]
        #print(" Edges out: %s" % edges_out)
        for next_node in edges_out:
            #print("   Going to node %d" % next_node)
            self._traverse_from_node(next_node, node, n_bp_traversed)


    def align(self):
        nodes = list(sorted(self._nodes))
        edges = sorted(self._edges, key=lambda e: e[0])
        sequences = [self._sequence_graph.get_sequence_on_directed_node(node) for node in nodes]
        if self.is_reverse:
            nodes = [-n for n in nodes]
            nodes = sorted(nodes)
            edges = [(-e[0], -e[1]) for e in edges]

        #print("  --- ")
        #print(self.sequence)
        #print(len(self.sequence))
        #print(nodes)
        #print(edges)
        #print(sequences)
        aligned_to_nodes, score = align(nodes, sequences, edges, self.sequence)
        if self.is_reverse:
            aligned_to_nodes = [-n for n in aligned_to_nodes]
        return aligned_to_nodes, score


def align(nodes, node_sequences, edges, sequence):
    # nodes is list of node ids
    # node_sequences is list of node sequences
    # edges is list of tuples (from node, to node)

    match = 1
    mismatch = 4
    gap_open = 6
    gap_extension = 1
    read_seq = "TTTTTTT"

    nttable = gssw_create_nt_table()
    mat = gssw_create_score_matrix(match, mismatch)
    gssw_nodes = {}
    for i, node in enumerate(nodes):
        n = gssw_node_create(None, node, node_sequences[i], nttable, mat)
        gssw_nodes[node] = n

    for edge in edges:
        gssw_nodes_add_edge(gssw_nodes[edge[0]], gssw_nodes[edge[1]])

    graph = gssw_graph_create(len(nodes))
    for node in gssw_nodes.values():
        gssw_graph_add_node(graph, node)

    gssw_graph_fill(graph, sequence, nttable, mat, gap_open, gap_extension, 0, 0, 15, 2, True)
    mapping = test_wrapper(graph,
                           sequence,
                           len(sequence),
                           nttable,
                           mat,
                           gap_open,
                           gap_extension,
                           0, 0)

    score = mapping[-1]
    return mapping[0:-1], score*2


if __name__ == "__main__":
    import sys

    """
    nodes = [534647, 534648, 534649, 534650, 534651, 534652, 534653, 534654, 534655, 534656, 534657, 534658, 534659, 534660, 534661, 534662, 534663, 534664, 534665, 534666, 534667]
    edges = [(534647, 534649), (534647, 534648), (534648, 534650), (534649, 534650), (534650, 534651), (534650, 534652),
             (534651, 534653), (534652, 534653), (534653, 534654), (534653, 534655), (534654, 534656), (534655, 534656),
             (534656, 534659), (534656, 534657), (534656, 534658), (534657, 534658), (534658, 534659), (534659, 534660),
             (534659, 534662), (534659, 534661), (534660, 534662), (534661, 534662), (534661, 534660), (534662, 534665),
             (534662, 534663), (534663, 534666), (534663, 534664), (534664, 534666), (534665, 534664), (534665, 534666),
             (534666, 534667)]
    sequences = ['gaa', 'a', 'g', 'aa', 'g', 'a', 'aagaagaagaagaggaggagg', 'a', 'g', 'ggagggagaa', 'aagg', 'aagg', 'aaggaaggaaggaagg', 'agg', 'a', 'aaaaa', 'a', 'c', 'c', 'cccacattgcagactcctatttgaggaagctg', 'acctctacaatc']
    alignment, score = align(nodes, sequences, edges, "CCCCC")
    print(alignment, score)
    """
    nodes, score = align(
        [1, 2, 3, 4], ["AGG", "AAAA", "A", "A"], [(1, 2), (1, 3), (3, 2), (3, 4), (2, 4)], "CCCCCC"
    )

    nodes, score = align(
        [59, 60, 61, 62], ["AGG", "AAAA", "A", "A"], [(59, 60), (59, 61), (61, 60), (61, 62), (60, 62)], "CCCCCC"
    )
    print(nodes, score)
    sys.exit()
    nodes, score = align([1, 2, 3, 4], ["G", "AAAAA", "C", "TTTT"], [(1, 2), (1, 3), (2, 4), (3, 4)], "GCTTT")
    sys.exit()
    for i in range(0, 100000):
        read = "TTATCATGGGATCTATATACTGATCTACCTTTGGTGACTGTACGTACGTAGCA"
        nodes = align([1, 2, 3, 4, 5, 6], ["AAAACCAGTTATCACG", "CCCCCCTTCTTCTCTTTCGACGACTA", "GGATCTATATACT", "GATCT", "GTGACTGTACGTACGTAGCA", "GGG"], [(1, 2), (1, 3), (2, 4), (3, 4), (4, 5), (5, 6), (4, 6)], read)
        #print(nodes)

    sys.exit()
    nodes, score = align(
        [1, 3, 2, 4], ["AGG", "AAAA", "A", "A"], [(1, 3), (3, 2), (1, 4), (2, 4), (1, 4), (3, 4)], "CCCCCC"
    )
    print(nodes, score)
    sys.exit()
    nodes, score = align([1, 2, 3, 4], ["G", "AAAAA", "C", "TTTT"], [(1, 2), (1, 3), (2, 4), (3, 4)], "GCTTT")
    sys.exit()
    for i in range(0, 100000):
        read = "TTATCATGGGATCTATATACTGATCTACCTTTGGTGACTGTACGTACGTAGCA"
        nodes = align([1, 2, 3, 4, 5, 6], ["AAAACCAGTTATCACG", "CCCCCCTTCTTCTCTTTCGACGACTA", "GGATCTATATACT", "GATCT", "GTGACTGTACGTACGTAGCA", "GGG"], [(1, 2), (1, 3), (2, 4), (3, 4), (4, 5), (5, 6), (4, 6)], read)
        #print(nodes)

