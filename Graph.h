/*
 * Graph.h
 *
 *  Created on: 5Dec.,2017
 *      Author: Lijun Chang
 *      Email: ljchang@outlook.com
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include "Defines.h"
#include "Utility.h"
#include "Timer.h"
#include "ListLinearHeap.h"
#include "UnionFind.h"

struct Edge {
	Edge *pre, *next;
	Edge *reverse;
	ui vertex;
};

struct Stack_node {
	ui vertex, parent_cid;
	ui L, H;
	char level;
};

class Graph {
private:
	std::string dir; //input graph directory
	ui n; //#nodes of the graph
	ui m; //#edges of the graph

	ui *pstart; //start positions of neighbors of vertices in the array "edges"
	ui *edges; //concatenation of neighbors of all vertices

public:
	Graph(const char *_dir) ;
	~Graph() ;

	void read_graph_binary() ;

	// compute k-edge connected components
	void k_edge_connected_component(ui K, std::string output_file) ;
	// compute k-edge connected components in a space-effective manner
	void k_edge_connected_component_space(ui K, std::string output_file) ;
	// edge-connectivity-based graph decomposition, divide and conquer
	void edge_connectivity_decomposition_DCs(bool mspt, std::string output_file, std::string eccsizes_file) ;
	// edge-connectivity-based graph decomposition, bottom-up (large to small value), core-optimized
	void edge_connectivity_decomposition_BUso(bool mspt, std::string output_file, std::string eccsizes_file) ;
	// edge-connectivity-based graph decomposition, bottom-up (large to small value)
	void edge_connectivity_decomposition_BUs(bool mspt, std::string output_file, std::string eccsizes_file) ;
	// edge-connectivity-based graph decomposition, top-down (small to large value)
	void edge_connectivity_decomposition_TDs(bool mspt, std::string output_file, std::string eccsizes_file) ;

private:
	void to_hierarchy_tree(std::vector<std::pair<std::pair<ui,ui>, ui> > vpp, std::string output_file, std::string eccsizes_file) ;
	// core decomposition without invoking data structures
	ui core_decomposition(ui *peel_sequence, ui *core) ;
	// k-core based pruning
	void k_core_prune(ui K, ui *Q, ui Q_n, char *computed, ui *degree, ui *pend) ;
	void print_kecc(std::string output_file, ui K, ui c_n, ui *cstart, ui *ids) ;

	void get_degrees(const ui K, const ui *active_ids, const ui active_n, ui *degree, UnionFind *UF, const ui *representative, const ui *pend_global, ui *pend, const ui *core, ui *adj_next_global, ui *adj_next) ;

	void contract_graph(ui pos, ui *adj_next, const std::vector<std::pair<ui,ui> > &contractions, ui *adj_last_local) ;
	ui kECC(const ui K, ui *pend, const ui *active_ids, const ui active_n, ui *cstart, ui *ids, ui *degree, ui *Q, char *vis, char *computed, ui *pend_local, UnionFind *UF, const ui *representative, UnionFind *UF_local, ui *representative_local, ui *adj_next, ui *adj_next_local, ui *adj_last_local, ui *sv_next_local, ui *sv_last_local, ui *keys, ListLinearHeap *heap) ;
	ui get_active_ids(const ui s, const ui *active_component, ui *active_ids, ui *degree, char *vis, ui *pend, const ui *pend_global, UnionFind *UF, const ui *representative, ui *adj_next, ui *adj_next_child) ;
	ui k_core_prune(const ui K, const ui *pend, ui *Q, ui Q_n, char *computed, ui *degree, UnionFind *UF, const ui *representative, const ui *adj_next) ;
	ui initialize_pgraph(const ui s, ui *pend, ui *Q, char *vis, char *computed, UnionFind *UF, const ui *representative, ui *adj_next, ui *pend_local, UnionFind *UF_local, ui *representative_local, ui *sv_next_local, ui *sv_last_local, ui *adj_next_local, ui *adj_last_local) ;
	ui remove_inter_edges(const ui K, const ui c_n, const ui new_cn, const ui *cstart, const ui *ids, ui *pend, ui *Q, char *computed, ui *degree, UnionFind *UF_local, const ui *representative_local, ui *adj_next) ;
	void print_edge_connectivities(const char *alg, std::vector<std::pair<std::pair<ui,ui>, ui> > &vpp, bool mspt) ;
	void print_active_graph(const ui *active_ids, const ui active_n, const ui *pend, const ui *adj_next, UnionFind *UF, const ui *representative) ;
	void print_keccs(ui level, ui K, ui c_n, ui *cstart, ui *ids) ;

	// initialize the partition graph
	void initialize_pgraph(ui s,ui *Q, char *vis, char *computed, ui *pend, ui *pend_local, UnionFind *UF, ui *representative, ui *sv_next, ui *sv_last, ui *adj_next, ui *adj_last) ;
	ui decomposition(const ui s, const ui K, ui *cstart, ui *ids, ui c_n, ui *Q, ui *keys, char *vis, ui *pend_local, ListLinearHeap *heap, ui *sv_next, ui *sv_last, ui *adj_next, ui *adj_last, UnionFind *UF, ui *representative) ;
	void remove_inter_edges(ui K, ui c_n, ui new_cn, ui *cstart, ui *ids, ui *pend, ui *Q, char *computed, ui *degree, UnionFind *UF) ;

	// construct the partition graph
	void construct_pgraph(ui s, ui *Q, char *vis, char *computed, ui *pend, Edge **graph_head, Edge *graph_edges, ui *sv_next, ui *sv_last) ;
	ui decomposition(ui s, ui K, ui *cstart, ui *ids, ui c_n, ui *Q, ui *keys, char *vis, Edge **graph_head, Edge *graph_edges, ListLinearHeap *heap, ui *sv_next, ui *sv_last) ;
	void merge(Edge **graph_head, Edge *edges, ui u, ui v, ui *keys, ui *sv_next, ui *sv_last) ;
	void add_edge(ui u, Edge *e, Edge **graph_head) ;
	void delete_edge(ui u, Edge *e, Edge **graph_head) ;
	void remove_inter_edges(ui K, ui c_n, ui new_cn, ui *cstart, ui *ids, ui *component_id, ui *pend, ui *Q, char *computed, ui *degree) ;
};

#endif /* GRAPH_H_ */
