import logging
from pygssw.gssw import gssw_node_create, gssw_create_nt_table, gssw_create_score_matrix, \
    gssw_graph_create, gssw_create_score_matrix, gssw_nodes_add_edge, gssw_graph_add_node, \
    gssw_graph_fill, gssw_graph_print_score_matrices, gssw_graph_trace_back, test_wrapper


class OffsetBasedGraphAligner:
    def __init__(self, ob_graph, sequence_graph, start_node, sequence,
                 n_bp_to_traverse_right=None, n_bp_to_traverse_left=None):
        self.ob_graph = ob_graph
        self.start_node = start_node
        self._sequence_graph = sequence_graph
        self.sequence = sequence

        if n_bp_to_traverse_right is None:
            assert n_bp_to_traverse_left is None, "Both right and left bp must be set if either is set"
            self._n_bp_to_traverse_right = len(sequence) + 10
            self._n_bp_to_traverse_left = 0
        else:
            self._n_bp_to_traverse_left = n_bp_to_traverse_left
            self._n_bp_to_traverse_right = n_bp_to_traverse_right

        self.is_reverse = False
        #self.adj_list = ob_graph.adj_list
        if self.start_node < 0:
            self.is_reverse = True
            #self.adj_list = ob_graph.reverse_adj_list


        self._nodes = set()
        self._edges = set()

        self._traversed_nodes = {}  # node id => bp to that node

        self._find_all_local_nodes_and_edges()

    def _find_all_local_nodes_and_edges(self):
        self._traverse_from_node(self.start_node, None, n_bp_traversed=0, is_traversing_left=False)
        self._traverse_from_node(-self.start_node, None, n_bp_traversed=0, is_traversing_left=True)

    def _traverse_from_node(self, node, prev_node, n_bp_traversed, is_traversing_left=False):
        adj_list = self.ob_graph.adj_list
        n_bp_to_traverse = self._n_bp_to_traverse_right
        if is_traversing_left:
            n_bp_to_traverse = self._n_bp_to_traverse_left
            adj_list = self.ob_graph.reverse_adj_list

        #print("On node %d. Prev: %s" % (node, prev_node))
        node = int(node)
        if prev_node is not None:
            #print("Adding edge %d-%d" % (prev_node, node))
            if (is_traversing_left and abs(prev_node) < abs(node)) or (not is_traversing_left and node < prev_node):
                #logging.warning("Graph is not sorted. Skipping edges that are not sorted when locally aligning")
                pass
            else:
                if is_traversing_left:
                    self._edges.add((abs(node), abs(prev_node)))
                else:
                    self._edges.add((prev_node, node))

        if node in self._traversed_nodes:  # and self._traversed_nodes[node] <= n_bp_traversed:
            #print("   Stopping. already processed")
            # Already processed this node with a possibly shorter path
            return

        self._traversed_nodes[node] = n_bp_traversed

        if is_traversing_left:
            self._nodes.add(abs(node))
        else:
            self._nodes.add(node)


        n_bp_traversed += self.ob_graph.blocks[node].length()

        if n_bp_traversed > n_bp_to_traverse:
            #print("    Stopping because length (traversed: %d, max: %d)" % (n_bp_traversed, n_bp_to_traverse))
            return

        edges_out = adj_list[node]
        #print(" Edges out: %s" % edges_out)
        for next_node in edges_out:
            #print("   Going to node %d" % next_node)
            self._traverse_from_node(next_node, node, n_bp_traversed, is_traversing_left=is_traversing_left)


    def align(self):
        nodes = list(sorted(self._nodes))
        edges = sorted(self._edges, key=lambda e: e[0])
        sequences = [self._sequence_graph.get_sequence_on_directed_node(node) for node in nodes]
        if self.is_reverse:
            # Hacky conversion, we need nodes to increase
            max_node = max(abs(n) for n in nodes)
            nodes = [max_node - abs(n) for n in nodes]
            nodes = sorted(nodes)
            edges = [(max_node - abs(e[0]), max_node - abs(e[1])) for e in edges]

        aligned_to_nodes, score = align(nodes, sequences, edges, self.sequence)
        if self.is_reverse:
            # Convert node ids back
            aligned_to_nodes = [-n + max_node for n in aligned_to_nodes]
        return aligned_to_nodes, score

def compress_node_ids(nodes, edges):

    node_id = 1
    forward_mapping = {}
    back_mapping = {}

    new_nodes = []
    new_edges = []
    for node in nodes:
        forward_mapping[node] = node_id
        back_mapping[node_id] = node
        new_nodes.append(node_id)
        node_id += 1

    for edge in edges:
        edge = (forward_mapping[edge[0]], forward_mapping[edge[1]])
        new_edges.append(edge)

    return new_nodes, new_edges, back_mapping


def decompress_node_ids(nodes, mapping):
    new_nodes = []
    for node in nodes:
        new_nodes.append(mapping[node])

    return new_nodes 


def align(nodes, node_sequences, edges, sequence):
    # nodes is list of node ids
    # node_sequences is list of node sequences
    # edges is list of tuples (from node, to node)


    #max_node = min(nodes) - 1
    #nodes = [n - max_node for n in nodes]
    #edges = [(e[0] - max_node, e[1] - max_node) for e in edges]
    nodes, edges, node_mapping = compress_node_ids(nodes, edges)
    #print(nodes)
    #print(edges)
    
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

    #gssw_graph_fill(graph, sequence, nttable, mat, gap_open, gap_extension, 0, 0, 0, 2, True)
    mapping = test_wrapper(graph,
                           sequence,
                           len(sequence),
                           nttable,
                           mat,
                           gap_open,
                           gap_extension,
                           0, 0)

    score = mapping[-1]
    aligned_nodes = mapping[0:-1]
    #aligned_nodes = [n + max_node for n in aligned_nodes]
    aligned_nodes = decompress_node_ids(aligned_nodes, node_mapping)
    return aligned_nodes, score*2


if __name__ == "__main__":
    import sys
    read = "AAA"
    nodes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 16]
    edges = [(1, 2), (2, 3), (3, 4), (3, 5), (4, 6), (5, 6), (6, 7), (7, 8), (8, 9), (9, 16), (9, 12)]
    sequences = ["A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A"]
    alignment, score = align(nodes, sequences, edges, read)
    print(alignment, score)
    sys.exit()

    read = "GCGGTTAGGAGCTGGAGACCAGCCTTGCCAACCCAGCGAAACCCCGTCCCCACCAAAAAAATACGAAAAGCAGTCAGGCGTGGCGGCGTGCGCCTGCAATCGCAGGCACTCGGCAGGTTGAGGCAGGAGAATCAGGCAGGGAGGTTGCAG"
    nodes = [1731517, 1731518, 1731519, 1731520, 1731521, 1731522, 1731523, 1731524, 1731525, 1731526, 1731527, 1731528,
     1731529, 1731530, 1731531, 1731532, 1731533, 1731534, 1731535]
    edges = [(1731517, 1731518), (1731518, 1731519), (1731519, 1731520), (1731520, 1731521), (1731521, 1731522),
     (1731522, 1731523), (1731522, 1731524), (1731523, 1731525), (1731524, 1731525), (1731525, 1731526),
     (1731525, 1731527), (1731526, 1731528), (1731527, 1731528), (1731528, 1731529), (1731529, 1731530),
     (1731529, 1731531), (1731530, 1731532), (1731531, 1731532), (1731532, 1731533), (1731533, 1731534),
     (1731534, 1731535)]
    sequences = ['gtgaacgagagtccgtctgcaatcccggcacc', 'tcgggaggccgaggctggcggatcactcgcgg', 'ttaggagctgcagaccagcccagccaacacag',
     'cgaatccccgtctccaccaaaaaaatacgaaa', 'accagtcaggcgtggcggcgtgcgcctgcaat', 'cgcaggca', 't', 'c', 'tcggcaa', 'a', 'g',
     'ctgaggcaagagaatcaggcagggaggttgca', 'gtgagccgagat', 'a', 'g', 'gcagcagtaccgtccagctttggctcggcatc',
     'agagggagaccgtggaaagagagggagtggga', 'gaccgtggggagagggagagggagagggaggg', 'ggagggggagggagagggagagggagagggag']
    alignment, score = align(nodes, sequences, edges, read)
    print(alignment, score)

    sys.exit()

    read = "AAA"
    nodes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    edges = [(1, 2), (2, 3), (3, 4), (3, 5), (4, 6), (5, 6), (6, 7), (7, 8), (8, 9), (9, 10), (9, 11)]
    sequences = ["A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A"]
    alignment, score = align(nodes, sequences, edges, read)
    print(alignment, score)
    sys.exit()


    nodes = [1731517, 1731518, 1731519, 1731520, 1731521, 1731522, 1731523, 1731524, 1731525, 1731526, 1731527, 1731528,
     1731529, 1731530, 1731531, 1731532, 1731533, 1731534, 1731535]
    edges = [(1731517, 1731518), (1731518, 1731519), (1731519, 1731520), (1731520, 1731521), (1731521, 1731522),
     (1731522, 1731523), (1731522, 1731524), (1731523, 1731525), (1731524, 1731525), (1731525, 1731526),
     (1731525, 1731527), (1731526, 1731528), (1731527, 1731528), (1731528, 1731529), (1731529, 1731530),
     (1731529, 1731531), (1731530, 1731532), (1731531, 1731532), (1731532, 1731533), (1731533, 1731534),
     (1731534, 1731535)]


    #nodes = [n - 1731516 for n in nodes]
    #edges = [(e[0] - 1731516, e[1] - 1731516) for e in edges]
    sequences = ['gtgaacgagagtccgtctgcaatcccggcacc', 'tcgggaggccgaggctggcggatcactcgcgg', 'ttaggagctgcagaccagcccagccaacacag',
     'cgaatccccgtctccaccaaaaaaatacgaaa', 'accagtcaggcgtggcggcgtgcgcctgcaat', 'cgcaggca', 't', 'c', 'tcggcaa', 'a', 'g',
     'ctgaggcaagagaatcaggcagggaggttgca', 'gtgagccgagat', 'a', 'g', 'gcagcagtaccgtccagctttggctcggcatc',
     'agagggagaccgtggaaagagagggagtggga', 'gaccgtggggagagggagagggagagggaggg', 'ggagggggagggagagggagagggagagggag']
    alignment, score = align(nodes, sequences, edges, "CCCCC")
    print(alignment, score)

    sys.exit()
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
        [1, 2, 3, 4], ["AGG", "AAAA", "A", "A"], [(1, 2), (1, 3), (3, 2), (3, 4), (2, 4)], read
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

