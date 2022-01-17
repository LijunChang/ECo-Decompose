# Edge Connectivity-based Graph Decomposition Algorithms
This repository implements an efficient algorithm for edge connectivity-based hierarchical graph decomposition.

# Graph Datasets
Each graph is presented by two binary files, b_adj.bin and b_degree.bin (e.g. datasets/CA-GrQc/b_adj.bin and datasets/CA-GrQc/b_degree.bin).

In the b_degree.bin, the first line is a single number checking whether the size of unsigned int in bytes of the machine is consistent with the binary files.

The second line is a single number representing the number of vertices (n) of the graph.

The third line is a single number representing the number of directed edges (2*m) of the graph (each undirected edge counts as two directed edges).

For the next n lines, each contains a single number corrsponding to the degree of a vertex (e.g. the next first line contains the degree of vertex 0).

In the b_adj.bin, there are n lines in total.

Each line contains multiple numbers representing the neighbours of a vertex (e.g. the first line includes neighbours of vertex 0). 

## Compile

```
make
```
It generates an executable "eco_decompose"

## 1. Run kecc-space

```
./eco_decompose -g datasets/CA-GrQc/ -a kecc-space -k 10
```

## 2. Run eco-dcs

```
./eco_decompose -g datasets/CA-GrQc/ -a eco-decompose-dcs
```

