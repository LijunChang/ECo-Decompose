This repository implements edge connectivity-based hierarchical graph decomposition algorithms proposed in our VLDB 2022 paper. If you are using the code, please cite our paper.
<pre>
Lijun Chang and Zhiyi Wang.
[pdf/vldb22.pdf](A Near-Optimal Approach to Edge Connectivity-Based Hierarchical Graph Decomposition.)
Proc. VLDB Endow. 15(6), (2022)
</pre>

## Compile the code

```sh
$ make clean
$ make
```
It generates an executable "eco_decompose".

## Run the code

Different algorithms can be invoked by executing "eco_decompose". You can find how to use the code by
```sh
$ ./eco_decompose -h
```

An example of computing edge connectivity-base hierarchical graph decomposition for the dataset CA-GrQc by the ECo-DC-AA algorithm is as follows
```sh
$ ./eco_decompose -g datasets/CA-GrQc/ -a eco-decompose-dcs -o result.txt
```

An example of computing 10-edge connected components for the dataset CA-GrQc by the algorithm KECC-AA is as follows
```sh
$ ./eco_decompose -g datasets/CA-GrQc/ -a kecc-space -k 10 -o result.txt
```

## Data format
Each graph is represented by two binary files, b_adj.bin and b_degree.bin (e.g. datasets/CA-GrQc/b_adj.bin and datasets/CA-GrQc/b_degree.bin). More details of the data format can be found in [https://github.com/LijunChang/Cohesive_subgraph_book/tree/master/datasets](https://github.com/LijunChang/Cohesive_subgraph_book/tree/master/datasets)


[//]: # "In the b_degree.bin, the first line is a single number checking whether the size of unsigned int in bytes of the machine is consistent with the binary files."

[//]: # "The second line is a single number representing the number of vertices (n) of the graph."

[//]: # "The third line is a single number representing the number of directed edges (2*m) of the graph (each undirected edge counts as two directed edges)."

[//]: # "For the next n lines, each contains a single number corrsponding to the degree of a vertex (e.g. the next first line contains the degree of vertex 0)."

[//]: # "In the b_adj.bin, there are n lines in total."

[//]: # "Each line contains multiple numbers representing the neighbours of a vertex (e.g. the first line includes neighbours of vertex 0)."

