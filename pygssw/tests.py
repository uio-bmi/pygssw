from align import Aligner
from offsetbasedgraph import Graph, Block

def test_get_nodes_and_edges():
    graph = Graph({1: Block(10), 2: Block(10), 3: Block(5), 4: Block(10), 5: Block(10)},
                   {1: [2, 3], 2: [4], 3: [4], 4: [5]})

    read =  "A" * 34
    aligner = Aligner(graph, None, 1, read)
    print(aligner._nodes)
    print(aligner._edges)


if __name__ == "__main__":
    test_get_nodes_and_edges()