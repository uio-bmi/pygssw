from gssw import gssw_node_create, gssw_create_nt_table, gssw_create_score_matrix, \
    gssw_graph_create, gssw_create_score_matrix, gssw_nodes_add_edge, gssw_graph_add_node, \
    gssw_graph_fill, gssw_graph_print_score_matrices, gssw_graph_trace_back, test_wrapper
import numpy as np

match = 1
mismatch = 4
gap_open = 6
gap_extension = 1
read_seq = "TTTTTTT"

nttable = gssw_create_nt_table()
mat = gssw_create_score_matrix(match, mismatch)
nodes = []
nodes.append(gssw_node_create(None, 10, "CCACT", nttable, mat))
nodes.append(gssw_node_create(None, 2, "TTT", nttable, mat))
nodes.append(gssw_node_create(None, 3, "TTTCC", nttable, mat))
gssw_nodes_add_edge(nodes[0], nodes[1])
gssw_nodes_add_edge(nodes[1], nodes[2])
graph = gssw_graph_create(3)
gssw_graph_add_node(graph, nodes[0])
gssw_graph_add_node(graph, nodes[1])
gssw_graph_add_node(graph, nodes[2])

gssw_graph_fill(graph, read_seq, nttable, mat, gap_open, gap_extension, 0, 0, 15, 2, True)
mapping = test_wrapper(graph,
                            read_seq,
                            len(read_seq),
                            nttable,
                            mat,
                            gap_open,
                            gap_extension,
                            0, 0)

print(mapping)
print("Done")
#gssw_graph_print_score_matrices(graph, "AAA", len("AAA"), None)
print("...")

"""
nodes[0] = (gssw_node*)gssw_node_create("A", 1, ref_seq_1, nt_table, mat);
    nodes[1] = (gssw_node*)gssw_node_create("B", 2, ref_seq_2, nt_table, mat);
    nodes[2] = (gssw_node*)gssw_node_create("C", 3, ref_seq_3, nt_table, mat);
    nodes[3] = (gssw_node*)gssw_node_create("D", 4, ref_seq_4, nt_table, mat);
    
    // makes a diamond
    gssw_nodes_add_edge(nodes[0], nodes[2]);
    gssw_nodes_add_edge(nodes[1], nodes[3]);
    gssw_nodes_add_edge(nodes[2], nodes[3]);
    
    gssw_graph* graph = gssw_graph_create(4);
    //memcpy((void*)graph->nodes, (void*)nodes, 4*sizeof(gssw_node*));
    //graph->size = 4;
    gssw_graph_add_node(graph, nodes[0]);
    gssw_graph_add_node(graph, nodes[1]);
    gssw_graph_add_node(graph, nodes[2]);
    gssw_graph_add_node(graph, nodes[3]);
    
    gssw_graph_fill(graph, read_seq, nt_table, mat, gap_open, gap_extension, 0, 0, 15, 2, true);
    gssw_graph_print_score_matrices(graph, read_seq, strlen(read_seq), stdout);
    gssw_graph_mapping* gm = gssw_graph_trace_back (graph,
                                                    read_seq,
                                                    strlen(read_seq),
                                                    nt_table,
                                                    mat
                                                    gap_open,
                                                    gap_extension,
                                                    0, 0);

    printf("Optimal local mapping:\n");
    gssw_print_graph_mapping(gm, stdout);
    gssw_graph_mapping_destroy(gm);
"""

print(nttable)