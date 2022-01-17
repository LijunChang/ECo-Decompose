# Edge Connectivity-based Graph Decomposition Algorithms
This repository implements an efficient algorithm for edge connectivity-based hierarchical graph decomposition.

# Graph Datasets
Each graph is presented by two binary files, b_adj.bin and b_degree.bin (e.g. datasets/CA-GrQc/b_adj.bin and datasets/CA-GrQc/b_degree.bin).

In the b_degree.bin, the first line is a single number checking whether the size of unsigned int in bytes of the machine is consistent with the binary files.

The second line is a single number representing the number of vertices (n) of the graph.

The third line is a single number representing the number of directed edges (2*m) of the graph (each undirected edge counts as two directed edges).

For the next n lines, each contains a single number corrsponding to the degree of a vertex (e.g. the next first line contains the degree of vertex 0).

In the b_adj.bin, there are n lines.

Each line contains multiple numbers representing the neighbours of a vertex (e.g. the first line includes neighbours of vertex 0). 

## Compile

```
make
```
It generates an executable "eco_decompose"

## 1. Run kecc or kecc-space

```
./eco_decompose ../datasets/as-skitter/ kecc 10 output
./eco_decompose ../datasets/as-skitter/ kecc-space 10 output
```
Note that, the fourth parameter is an integer that specifies the value of k. kecc-space is more space effective than kecc; that is, kecc-space consumes smaller main memory space.

This implements the algorithm proposed in the following SIGMOD'13 paper, which computes all k-edge connected components of a graph for a given k.

Lijun Chang, Jeffrey Xu Yu, Lu Qin, Xuemin Lin, Chengfei Liu, and Weifa Liang <br/>
**Efficiently Computing k-Edge Connected Components via Graph Decomposition** <br/>
*Proceedings of the ACM SIGMOD International Conference on Management of Data* (SIGMOD’13), 2013

