/*
 * Graph.cpp
 *
 *  Created on: 5Dec.,2017
 *      Author: Lijun Chang
 *      Email: ljchang@outlook.com
 */

#include "Graph.h"

using namespace std;

Graph::Graph(const char *_dir) {
	dir = string(_dir);

	n = 0;
	m = 0;

	pstart = nullptr;
	edges = nullptr;
}

Graph::~Graph() {
	if(pstart != nullptr) {
		delete[] pstart;
		pstart = nullptr;
	}
	if(edges != nullptr) {
		delete[] edges;
		edges = nullptr;
	}
}

void Graph::edge_connectivity_decomposition_DCs(bool mspt, string output_file, string eccsizes_file) {
//#ifdef NDEBUG
//	printf("*** eco_decomposition_dcs (Release): %s ***\n", dir.c_str());
//#else
//	printf("*** eco_decomposition_dcs (Debug): %s ***\n", dir.c_str());
//#endif
	if(!mspt) printf("!!! For memory consideration, mspt should be set to true\n");

	ui *Q = new ui[n];
	ui *active_component = new ui[n];
	ui max_core = core_decomposition(Q, active_component);
	printf("*\tmax_core: %u\n", max_core);

	char *vis = new char[n];
	memset(vis, 0, sizeof(char)*n);

	stack<Stack_node> working_stack;

	for(ui i = 0;i < n;i ++) if(!vis[i]) {
		ui Q_n = 1;
		Q[0] = i; vis[i] = 1;
		for(ui j = 0;j < Q_n;j ++) {
			ui u = Q[j];
			for(ui k = pstart[u];k < pstart[u+1];k ++) if(!vis[edges[k]]) {
				vis[edges[k]] = 1;
				Q[Q_n ++] = edges[k];
			}
		}

		for(ui j = 0;j < Q_n;j ++) active_component[Q[j]] = i;
		Stack_node node;
		node.level = 0;
		node.vertex = i;
		node.parent_cid = i;
		node.L = 2;
		node.H = max_core;
		working_stack.push(node);
	}
	memset(vis, 0, sizeof(char)*n);

	const int max_level = 30;
	ui *adj_next_buf[max_level];
	for(ui i = 0;i < max_level;i ++) adj_next_buf[i] = nullptr;
	adj_next_buf[0] = new ui[n];

	for(ui i = 0;i < n;i ++) adj_next_buf[0][i] = i;

	UnionFind *UF = new UnionFind(n);
	UF->init(n);
	ui *representative = new ui[n]; // representative vertex of super vertex
	for(ui i = 0;i < n;i ++) representative[i] = i;

	UnionFind *UF_local = new UnionFind(n);
	UF_local->init(0);
	ui *representative_local = new ui[n]; // representative vertex of super vertex
	ui *sv_next_local = new ui[n]; // next of super vertex
	ui *sv_last_local = new ui[n]; // last of super vertex
	ui *adj_next_local = new ui[n];
	ui *adj_last_local = new ui[n];

	ui *ids = new ui[n];
	ui *cstart = new ui[n+1];

	ui *active_ids = new ui[n];
	ui *pend = new ui[n];
	ui *pend_global = new ui[n];
	ui *pend_local = new ui[n];

	ui *degree = new ui[n];
	ui *keys = new ui[n];
	ListLinearHeap *heap = new ListLinearHeap(n, max_core);

	for(ui i = 0;i < n;i ++) pend_global[i] = pstart[i+1];

	char *computed = new char[n];

	vector<pair<pair<ui,ui>, ui> > vpp;
	UnionFind *UF_spt = nullptr;

	if(mspt) {
		vpp.reserve(n);
		UF_spt = new UnionFind(n);
		UF_spt->init(n);
	}
	else vpp.reserve(m/2);

	vector<pair<ui,ui> > contractions;
	ui current_pos[max_level];

	ui max_k = 0;

	while(!working_stack.empty()) {
		ui s = representative[UF->UF_find(working_stack.top().vertex)], parent_cid = working_stack.top().parent_cid;
		ui L = working_stack.top().L, H = working_stack.top().H;
		int level = working_stack.top().level;
		ui cid = active_component[s];

#ifndef NDEBUG
		//printf("*********************************************************\n");
		//printf("processing s: %u, cid: %u, parent_cid: %u, L: %u, H: %u, level: %d\n", s, cid, parent_cid, L, H, level);
#endif

		if(level >= max_level) {
			active_component[s] = parent_cid;
			working_stack.pop();
			continue;
		}

		ui M = (L+H+1)/2;
		assert(L <= M <= H);

		if(level <= 0) {
			level = -level;
			contract_graph(current_pos[level], adj_next_buf[level], contractions, adj_last_local);
			current_pos[level] = contractions.size();
		}

		assert(level >= 0&&level+1 < max_level);
		if(adj_next_buf[level+1] == nullptr) adj_next_buf[level+1] = new ui[n];

		// set adj_next_buf[level] and pend
		ui active_n = get_active_ids(s, active_component, active_ids, degree, vis, pend, pend_global, UF, representative, adj_next_buf[level], adj_next_buf[level+1]);

#ifndef NDEBUG
		//print_active_graph(active_ids, active_n, pend, adj_next_buf[level], UF, representative);
#endif
		// only change adj_next_buf[level+1], pend_local, pend, but not adj_next_buf[level]
		ui c_n = kECC(M, pend, active_ids, active_n, cstart, ids, degree, Q, vis, computed, pend_local, UF, representative, UF_local, representative_local, adj_next_buf[level+1], adj_next_local, adj_last_local, sv_next_local, sv_last_local, keys, heap);
		assert(cstart[c_n] == active_n);

#ifndef NDEBUG
		//printf("\n");
		//print_keccs(1, M, c_n, cstart, ids);
#endif

#ifndef NDEBUG
		//test the result of kECC
		for(ui i = 0;i < c_n;i ++) {
			assert(representative[UF->UF_find(ids[cstart[i]])] == ids[cstart[i]]);
			for(ui j = cstart[i];j < cstart[i+1];j ++) {
				assert(active_component[ids[j]] == cid);
				//printf("UF_local[%u]: %u, representative: %u\n", ids[j], UF_local->UF_find(ids[j]), representative_local[UF_local->UF_find(ids[j])]);
				assert(representative_local[UF_local->UF_find(ids[j])] == ids[cstart[i]]);
			}
		}
		for(ui i = 0;i < cstart[c_n];i ++) for(ui j = i+1;j < cstart[c_n];j ++) assert(ids[i] != ids[j]);

		//test the unmodificatiaon of the subgraph
		for(ui i = 0;i < active_n;i ++) {
			ui u = active_ids[i];
			ui tu = u;
			while(true) {
				assert(representative[UF->UF_find(tu)] == u);
				for(ui j = pstart[tu];j < pend[tu];j ++) assert(representative[UF->UF_find(edges[j])] != u);

				if(adj_next_buf[level][tu] == tu) break;
				else tu = adj_next_buf[level][tu];
			}
		}
#endif

#ifndef NDEBUG
		for(ui j = cstart[0];j < cstart[c_n];j ++) {
			ui u = ids[j];
			ui tu = u, ru = UF_local->UF_find(u);
			while(true) {
				for(ui k = pstart[tu];k < pend[tu];k ++) assert(UF_local->UF_find(edges[k]) == ru);
				if(adj_next_buf[level+1][tu] == tu) break;
				tu = adj_next_buf[level+1][tu];
			}
		}
#endif

		//printf("%d %d %d\n", L, M, H);
		if(M == H) {
			if(M > max_k) max_k = M;
			// use pend and adj_next_buf[level+1] to get the edges
			ui *adj_next = adj_next_buf[level+1];

			// assign M to edges
			for(ui j = cstart[0];j < cstart[c_n];j ++) {
				ui u = ids[j];
				ui tu = u;
				while(true) {
					for(ui k = pstart[tu];k < pend[tu];k ++) if(edges[k] > tu) {
						if(!mspt) {
							vpp.pb(make_pair(make_pair(tu, edges[k]), M));
#ifndef NDEBUG
							//printf("computed %u %u %u\n", tu, edges[k], M);
#endif
						}
						else if(UF_spt->UF_find(tu) != UF_spt->UF_find(edges[k])) {
							UF_spt->UF_union(tu, edges[k]);
							vpp.pb(make_pair(make_pair(tu, edges[k]), M));
						}
					}
					pstart[tu] = pend[tu];

					if(adj_next[tu] == tu) break;
					tu = adj_next[tu];
				}
			}

			// contract each M-edge connected component into a supervertex
			for(ui i = 0;i < c_n;i ++) {
				ui u = ids[cstart[i]];
				representative[UF->UF_find(u)] = u;
				for(ui j = cstart[i]+1;j < cstart[i+1];j ++) {
					representative[UF->UF_union(u, ids[j])] = u;
					contractions.pb(make_pair(u, ids[j]));
					//printf("inserted (%u,%u) to contractions\n", u, ids[j]);
				}
			}

#ifndef NDEBUG
			for(ui i = 0;i < c_n;i ++) {
				ui u = ids[cstart[i]];
				for(ui j = 0;j < n;j ++) if(representative[UF->UF_find(j)] == u) {
					for(ui k = pstart[j];k < pend_global[j];k ++) {
						//if(representative[UF->UF_find(edges[k])] == u) {
						//	printf("conflict %u, %u\n", j, edges[k]);
						//}
						assert(representative[UF->UF_find(edges[k])] != u);
					}
				}
			}
#endif

			if(L == M) {
				// assign M-1 to edges
				assert(cid == active_component[representative[UF->UF_find(s)]]);
				adj_next = adj_next_buf[level];
				for(ui i = 0;i < active_n;i ++) {
					ui u = active_ids[i];
					ui tu = u, ru = UF->UF_find(u);
					while(true) {
						for(ui &start = pstart[tu];start < pend_global[tu];start ++) {
							assert(UF->UF_find(edges[start]) != ru);
							if(active_component[representative[UF->UF_find(edges[start])]] != cid) {
#ifndef NDEBUG
								for(ui j = start;j < pend_global[tu];j ++) {
									assert(active_component[representative[UF->UF_find(edges[j])]] != cid);
								}
#endif
								break;
							}
							if(edges[start] > tu) {
								if(!mspt) {
									vpp.pb(make_pair(make_pair(tu, edges[start]), M-1));
#ifndef NDEBUG
									//printf("computed %u %u %u\n", tu, edges[start], M-1);
#endif
								}
								else if(UF_spt->UF_find(tu) != UF_spt->UF_find(edges[start])) {
									UF_spt->UF_union(tu, edges[start]);
									vpp.pb(make_pair(make_pair(tu, edges[start]), M-1));
								}
							}
						}

						if(adj_next[tu] == tu) break;
						tu = adj_next[tu];
					}
				}

				// contract all M-edge connected components into a single supervertex
				working_stack.pop();
				ui u = ids[0]; assert(cstart[0] == 0);
				representative[UF->UF_find(u)] = u;
				for(ui i = 1;i < c_n;i ++) {
					ui v = ids[cstart[i]];
					representative[UF->UF_union(u, v)] = u;
					contractions.pb(make_pair(u, v));
					//printf("inserted (%u,%u) to contractions\n", u, v);
				}

#ifndef NDEBUG
				for(ui i = 0;i < 1;i ++) {
					ui u = ids[cstart[i]];
					for(ui j = 0;j < n;j ++) if(representative[UF->UF_find(j)] == u) {
						for(ui k = pstart[j];k < pend_global[j];k ++) {
							assert(active_component[representative[UF->UF_find(edges[k])]] != cid);
							assert(representative[UF->UF_find(edges[k])] != u);
						}
					}
				}
#endif
				active_component[u] = parent_cid;
			}
			else {
				working_stack.top().H = M-1;
				level = working_stack.top().level;
				if(level > 0) working_stack.top().level = - level;
				else level = -level;
			}
		}
		else {
			working_stack.top().H = M-1;
			if(working_stack.top().level > 0) working_stack.top().level = - working_stack.top().level;
			level = working_stack.top().level;
			if(level > 0) working_stack.top().level = - level;
			else level = -level;

			ui pos = c_n;
			for(ui i = 0;i < c_n;i ++) {
				ui u = ids[cstart[i]];
				assert(representative[UF->UF_find(u)] == u);
				for(ui j = cstart[i];j < cstart[i+1];j ++) {
					assert(representative[UF->UF_find(ids[j])] == ids[j]);
					active_component[ids[j]] = u;
				}

				if(active_component[u] == cid) {
					assert(pos == c_n);
					pos = i;
					continue;
				}

				Stack_node node;
				node.vertex = u; node.parent_cid = cid; node.L = M+1; node.H = H; node.level = level+1;
				if(cstart[i+1] == cstart[i]+1) node.level = max_level;
				working_stack.push(node);
				//printf("\tpush %d %d\n", M+1, H);
			}

			if(pos != c_n) {
				Stack_node node;
				node.vertex = ids[cstart[pos]]; node.parent_cid = cid; node.L = M+1; node.H = H; node.level = level+1;
				if(cstart[pos+1] == cstart[pos]+1) node.level = max_level;
				//printf("\tpush %d %d\n", M+1, H);
				working_stack.push(node);
			}

			current_pos[level+1] = contractions.size();
		}
	}
	printf("*\tmax_k: %u\n", max_k);

	if(UF_spt != nullptr) delete UF_spt;

	for(ui i = 0;i < max_level;i ++) if(adj_next_buf[i] != nullptr) {
		delete[] adj_next_buf[i];
		adj_next_buf[i] = nullptr;
	}

	delete[] computed;

	delete[] degree;
	delete[] keys;
	delete heap;

	delete[] pend;
	delete[] pend_global;
	delete[] pend_local;

	delete[] active_ids;

	delete[] adj_next_local;
	delete[] adj_last_local;
	delete[] sv_next_local;
	delete[] sv_last_local;
	delete UF_local;
	delete[] representative_local;
	delete[] ids;
	delete[] cstart;

	delete UF;
	delete[] representative;

	delete[] vis;
	delete[] Q;
	delete[] active_component;

	to_hierarchy_tree(vpp, output_file, eccsizes_file);

	// print_edge_connectivities("eco-dcs", vpp, mspt);
}

void Graph::k_edge_connected_component_space(ui K, string output_file) {
//#ifdef NDEBUG
//	printf("*** k_edge_connected_component_space (Release): %s, %u ***\n", dir.c_str(), K);
//#else
//	printf("*** k_edge_connected_component_space (Debug): %s, %u ***\n", dir.c_str(), K);
//#endif

	if(K < 2) {
		printf("K must be at least 2!\n");
		return ;
	}

	ui *pend = new ui[n];
	for(ui i = 0;i < n;i ++) pend[i] = pstart[i+1];

	char *computed = new char[n];
	ui *Q = new ui[n], Q_n = 0;
	ui *degree = new ui[n];
	// k-core-based pruning
	memset(computed, 0, sizeof(char)*n);
	for(ui i = 0;i < n;i ++) {
		degree[i] = pend[i] - pstart[i];
		if(degree[i] < K) {
			Q[Q_n ++] = i;
			computed[i] = 1;
		}
	}
	k_core_prune(K, Q, Q_n, computed, degree, pend);

	UnionFind *UF = new UnionFind(n);
	ui *representative = new ui[n];
	ui *sv_next = new ui[n];
	ui *sv_last = new ui[n];
	ui *adj_next = new ui[n];
	ui *adj_last = new ui[n];

	ui *ids = new ui[n];
	ui *cstart = new ui[n+1];
	ui c_n = 0;
	cstart[0] = 0;

	ListLinearHeap *heap = new ListLinearHeap(n, K);
	ui *keys = new ui[n];

	char *vis = new char[n];
	memset(vis, 0, sizeof(char)*n);

	ui *pend_local = new ui[n];

	for(ui i = 0;i < n;) {
		if(computed[i]) {
			++ i;
			continue;
		}

		initialize_pgraph(i, Q, vis, computed, pend, pend_local, UF, representative, sv_next, sv_last, adj_next, adj_last);
		ui new_cn = decomposition(i, K, cstart, ids, c_n, Q, keys, vis, pend_local, heap, sv_next, sv_last, adj_next, adj_last, UF, representative);
		if(new_cn == c_n + 1) {
			for(ui j = cstart[c_n];j < cstart[c_n+1];j ++) computed[ids[j]] = 1;
			++ c_n;
		}
		else remove_inter_edges(K, c_n, new_cn, cstart, ids, pend, Q, computed, degree, UF);
	}

	//ui max_cc = 0;
	//for(ui i = 0;i < c_n;i ++) if(cstart[i+1] - cstart[i] > max_cc) max_cc = cstart[i+1]-cstart[i];
	//printf("maximum 2ecc size: %s\n", Utility::integer_to_string(max_cc).c_str());

	if(!output_file.empty()) print_kecc(output_file, K, c_n, cstart, ids);

	delete[] pend_local;
	delete[] vis;
	delete[] keys;
	delete[] pend;

	delete[] ids;
	delete[] cstart;

	delete[] representative;
	delete UF;
	delete[] sv_next;
	delete[] sv_last;
	delete[] adj_next;
	delete[] adj_last;

	delete[] computed;
	delete[] degree;
	delete[] Q;
}

// private member functions

void Graph::to_hierarchy_tree(vector<pair<pair<ui,ui>, ui> > vpp, string output_file, string eccsizes_file) {
	UnionFind *UF = new UnionFind(n);
	UF->init(n);

	int *vertex_to_eccnode = new int[n];
	for(ui i = 0;i < n;i ++) vertex_to_eccnode[i] = -1;

	ui ecc_n = 0;
	ui *parent = new ui[n];
	ui *weight = new ui[n];
	ui *representative_to_root = new ui[n];

	for(ui i = 0;i < vpp.size();i ++) {
		ui u = vpp[i].first.first, v = vpp[i].first.second, sc = vpp[i].second;
		if(vertex_to_eccnode[u] == -1&&vertex_to_eccnode[v] == -1) {
			parent[ecc_n] = ecc_n; representative_to_root[ecc_n] = ecc_n;
			vertex_to_eccnode[u] = vertex_to_eccnode[v] = ecc_n;
			weight[ecc_n] = sc;
			++ ecc_n;
		}
		else if(vertex_to_eccnode[u] == -1||vertex_to_eccnode[v] == -1) {
			if(vertex_to_eccnode[u] == -1) swap(u,v);
			ui ru = UF->UF_find(vertex_to_eccnode[u]);
			ui rru = representative_to_root[ru];
			if(weight[rru] != sc) {
				if(weight[rru] < sc) printf("WA in to_hierarchy_tree\n");
				parent[rru] = parent[ecc_n] = ecc_n;
				vertex_to_eccnode[v] = ecc_n;
				weight[ecc_n] = sc;
				ui new_r = UF->UF_union(ru, ecc_n);
				representative_to_root[new_r] = ecc_n;
				++ ecc_n;
			}
			else vertex_to_eccnode[v] = rru;
		}
		else {
			ui ru = UF->UF_find(vertex_to_eccnode[u]);
			ui rv = UF->UF_find(vertex_to_eccnode[v]);
			ui rru = representative_to_root[ru];
			ui rrv = representative_to_root[rv];
			if(ru == rv) {
				if(weight[rru] != sc) printf("!!! Hit\n");
				continue;
			}
			if(weight[rru] < sc||weight[rrv] < sc) printf("WA in to_hierarchy_tree\n");
			if(weight[rru] == sc||weight[rrv] == sc) {
				//if(weight[rru] == sc&&weight[rrv] == sc) printf("!!! Hit\n");
				if(weight[rrv] == sc) { swap(u,v); swap(ru, rv); swap(rru, rrv); }
				parent[rrv] = rru;
				ui new_r = UF->UF_union(ru,rv);
				representative_to_root[new_r] = rru;
			}
			else {
				parent[rru] = parent[rrv] = parent[ecc_n] = ecc_n;
				weight[ecc_n] = sc;
				ui new_r = UF->UF_union(ru, rv);
				new_r = UF->UF_union(new_r, ecc_n);
				representative_to_root[new_r] = ecc_n;
				++ ecc_n;
			}
		}
	}
	if(ecc_n > n) printf("WA in to_hierarchy_tree\n");

	bool bipartite = false;
	ui row_n = n;
	FILE *fbipartite = fopen((dir + string("/bipartite.txt")).c_str(), "r");
	if(fbipartite != NULL) {
		printf("* This is a bipartite graph!\n");
		bipartite = true;
		fscanf(fbipartite, "%u", &row_n);
		fclose(fbipartite);
	}

	ui *stack = new ui[n];
	ui stack_n = 0;
	char *in_stack = new char[n];
	memset(in_stack, 0, sizeof(char)*n);
	int *first_vertex = new int[n];
	for(ui i = 0;i < n;i ++) first_vertex[i] = -1;
	ui *ecc_size = representative_to_root;
	memset(ecc_size, 0, sizeof(ui)*n);
	ui *ecc_size_column = new ui[n];
	memset(ecc_size_column, 0, sizeof(ui)*n);
	for(ui i = 0;i < n;i ++) {
		if(vertex_to_eccnode[i] == -1) {
			printf("!!! There are isolated vertices in the input graph, please remove them\n");
			continue;
		}
		ui r = vertex_to_eccnode[i];
		ui end = n;
		if(first_vertex[r] != -1) {
			if(i < row_n) ++ ecc_size[first_vertex[r]];
			else ++ ecc_size_column[first_vertex[r]];
		}
		while(first_vertex[r] == -1) {
			ui tr = r;
			while(parent[r] != r&&weight[parent[r]] == weight[r]) r = parent[r];
			first_vertex[tr] = r;

			if(end == n) {
				if(i < row_n) ++ ecc_size[r];
				else ++ ecc_size_column[r];
			}

			if(in_stack[r]) break;
			in_stack[r] = 1;
			stack[-- end] = r;

			if(parent[r] == r) break;
			r = parent[r];
		}
		while(end < n) stack[stack_n ++] = stack[end ++];
	}
	for(ui i = stack_n;i > 0;i --) {
		ui u = stack[i-1];
		ui r = parent[u];
		if(r == u) continue;
		ecc_size[first_vertex[r]] += ecc_size[u];
		ecc_size_column[first_vertex[r]] += ecc_size_column[u];
	}

	if(!bipartite) {
		for(ui i = 0;i < n;i ++) if(ecc_size_column[i] > 0) {
			printf("WA for bipartite\n");
			break;
		}
	}

	if(!eccsizes_file.empty()&&!bipartite) {
		vector<pair<ui,ui> > vp;
		vp.reserve(stack_n);
		for(ui i = 0;i < stack_n;i ++) {
			ui pw = 0;
			if(parent[stack[i]] != stack[i]) pw = weight[parent[stack[i]]];
			for(ui j = pw+1;j <= weight[stack[i]];j ++) vp.pb(mp(j, ecc_size[stack[i]]));
		}
		sort(vp.begin(), vp.end());
		printf("*\tWriting ecc_sizes into file: %s\n", eccsizes_file.c_str());
		FILE *fout = Utility::open_file(eccsizes_file.c_str(), "w");
		for(ui i = 0;i < vp.size();i ++) fprintf(fout, "%u %u\n", vp[i].first, vp[i].second);
		fclose(fout);
	}
	if(!eccsizes_file.empty()&&bipartite) {
		vector<pair<ui,long long> > vp;
		vp.reserve(stack_n);
		for(ui i = 0;i < stack_n;i ++) {
			ui pw = 0;
			if(parent[stack[i]] != stack[i]) pw = weight[parent[stack[i]]];
			for(ui j = pw+1;j <= weight[stack[i]];j ++) vp.pb(mp(j, ((long long)ecc_size[stack[i]])*ecc_size_column[stack[i]]));
		}
		sort(vp.begin(), vp.end());
		printf("*\tWriting ecc_sizes into file: %s\n", eccsizes_file.c_str());
		FILE *fout = Utility::open_file(eccsizes_file.c_str(), "w");
		long long size = 0;
		for(ui i = 0;i < vp.size();i ++) {
			size += vp[i].second;
			if(i+1 == vp.size()||vp[i].first != vp[i+1].first) {
				fprintf(fout, "%u %lld\n", vp[i].first, size);
				size = 0;
			}
		}
		fclose(fout);
	}

	if(!output_file.empty()) {
		ui *id_map = representative_to_root;
		for(ui i = 0;i < n;i ++) id_map[i] = -1;
		for(ui i = 0;i < stack_n;i ++) id_map[stack[i]] = i;
		printf("*\tWriting hierarchy_tree into file: %s\n", output_file.c_str());
		FILE *fout = Utility::open_file(output_file.c_str(), "w");
		fprintf(fout, "%u %u\n", n, stack_n);
		for(ui i = 0;i < n;i ++) {
			if(first_vertex[vertex_to_eccnode[i]] == -1) fprintf(fout, "-1\n");
			else fprintf(fout, "%d\n", id_map[first_vertex[vertex_to_eccnode[i]]]);
		}
		for(ui i = 0;i < stack_n;i ++) {
			ui u = stack[i];
			int p = parent[u];
			if(p == u) p = -1;
			else p = id_map[first_vertex[p]];
			fprintf(fout, "%u %u %d\n", id_map[u], weight[u], p);
		}

		fclose(fout);
	}

	delete[] ecc_size_column;
	delete[] stack;
	delete[] first_vertex;
	delete[] in_stack;
	delete UF;
	delete[] parent;
	delete[] vertex_to_eccnode;
	delete[] weight;
	delete[] representative_to_root;
}

void Graph::get_degrees(ui K, const ui *active_ids, ui active_n, ui *degree, UnionFind *UF, const ui *representative, const ui *pend_global, ui *pend, const ui *core, ui *adj_next_global, ui *adj_next) {
	for(ui i = 0;i < active_n;i ++) {
		ui u = active_ids[i];
		ui tu = u, pre = u, pre_global = u;
		degree[u] = 0;
		while(true) {
			//printf("UF[%u]: %u, rep: %u, u: %u\n", tu, UF->UF_find(tu), representative[UF->UF_find(tu)], u);
			assert(representative[UF->UF_find(tu)] == u);
			if(pstart[tu] < pend_global[tu]) {
				ui &pend_t = pend[tu];
				pend_t = pstart[tu];
				while(pend_t < pend_global[tu]) {
					assert(representative[UF->UF_find(edges[pend_t])] != u);
					if(core != NULL&&core[edges[pend_t]] < K) break;
					else ++ pend_t;
				}
				if(pend_t > pstart[tu]) {
					degree[u] += pend_t-pstart[tu];
					adj_next[pre] = tu;
					pre = tu;
				}
			}
			ui next = adj_next_global[tu];
			adj_next_global[pre_global] = tu;
			pre_global = tu;
			if(next == tu) break;
			else tu = next;
		}
		adj_next[pre] = pre;
		adj_next_global[pre_global] = pre_global;
	}
}

void Graph::contract_graph(const ui pos, ui *adj_next, const vector<pair<ui,ui> > &contractions, ui *adj_last_local) {
#ifndef NDEBUG
	//printf("contract pairs");
	//for(ui i = pos;i < contractions.size();i ++) printf(" (%u,%u)", contractions[i].first, contractions[i].second);
	//printf("\n");
#endif

	for(ui i = pos;i < contractions.size();i ++) adj_last_local[contractions[i].first] = adj_last_local[contractions[i].second] = n;

	for(ui i = pos;i < contractions.size();i ++) {
		ui u = contractions[i].first;
		if(adj_last_local[u] == n) {
			ui tu = u;
			while(adj_next[tu] != tu) tu = adj_next[tu];
			adj_last_local[u] = tu;
		}
		u = contractions[i].second;
		if(adj_last_local[u] == n) {
			ui tu = u;
			while(adj_next[tu] != tu) tu = adj_next[tu];
			adj_last_local[u] = tu;
		}
	}

	for(ui i = pos;i < contractions.size();i ++) {
		ui u = contractions[i].first, v = contractions[i].second;
		adj_next[adj_last_local[u]] = v;
		adj_last_local[u] = adj_last_local[v];
	}
}

ui Graph::kECC(const ui K, ui *pend, const ui *active_ids, const ui active_n, ui *cstart, ui *ids, ui *degree, ui *Q, char *vis, char *computed, ui *pend_local, UnionFind *UF, const ui *representative, UnionFind *UF_local, ui *representative_local, ui *adj_next, ui *adj_next_local, ui *adj_last_local, ui *sv_next_local, ui *sv_last_local, ui *keys, ListLinearHeap *heap) {	ui Q_n = 0, c_n = 0;
	cstart[0] = 0;
	for(ui i = 0;i < active_n;i ++) {
		ui u = active_ids[i];
		computed[u] = 0;

		if(degree[u] < K) {
			Q[Q_n ++] = u;
			computed[u] = 1;
		}
	}
	Q_n = k_core_prune(K, pend, Q, Q_n, computed, degree, UF, representative, adj_next);
	for(ui i = 0;i < Q_n;i ++) {
		ids[cstart[c_n]] = Q[i];
		cstart[c_n+1] = cstart[c_n] + 1;
		++ c_n;
	}

	for(ui i = 0;i < active_n;) {
		ui u = active_ids[i];
		assert(representative[UF->UF_find(u)] == u);
		if(computed[u]) {
			++ i;
			continue;
		}

		ui size = initialize_pgraph(u, pend, Q, vis, computed, UF, representative, adj_next, pend_local, UF_local, representative_local, sv_next_local, sv_last_local, adj_next_local, adj_last_local);
#ifndef NDEBUG
		vector<ui> tmp;
		sort(Q, Q+size);
		for(ui j = 0;j < size;j ++) tmp.pb(Q[j]);
		for(ui j = 1;j < size;j ++) assert(tmp[j] > tmp[j-1]);
#endif
		ui new_cn = decomposition(u, K, cstart, ids, c_n, Q, keys, vis, pend_local, heap, sv_next_local, sv_last_local, adj_next_local, adj_last_local, UF_local, representative_local);

#ifndef NDEBUG
		print_keccs(2, K, new_cn-c_n, cstart+c_n, ids);
#endif

#ifndef NDEBUG
		for(ui j = cstart[c_n];j < cstart[new_cn];j ++) {
			char find = 0;
			for(ui k = 0;k < tmp.size();k ++) if(ids[j] == tmp[k]) find = 1;
			if(!find) printf("not found %u\n", ids[j]);
			assert(find);
		}
		tmp.clear();
		for(ui j = cstart[c_n];j < cstart[new_cn];j ++) tmp.pb(ids[j]);
		sort(tmp.begin(), tmp.end());
		for(ui j = 1;j < tmp.size();j ++) if(tmp[j] == tmp[j-1]) {
			printf("duplicate: %u\n", tmp[j]);
			assert(1 == 2);
		}
#endif

		assert(cstart[new_cn] == size + cstart[c_n]);

		if(new_cn == c_n + 1) {
			for(ui j = cstart[c_n];j < cstart[c_n+1];j ++) computed[ids[j]] = 1;
			++ c_n;
		}
		else {
			Q_n = remove_inter_edges(K, c_n, new_cn, cstart, ids, pend, Q, computed, degree, UF_local, representative_local, adj_next);
			Q_n = k_core_prune(K, pend, Q, Q_n, computed, degree, UF, representative, adj_next);
			for(ui j = 0;j < Q_n;j ++) {
				ids[cstart[c_n]] = Q[j];
				cstart[c_n+1] = cstart[c_n] + 1;
				++ c_n;
			}
		}
	}

	for(ui i = 0;i < c_n;i ++) if(cstart[i+1] == cstart[i] + 1) {
		ui u = ids[cstart[i]];
		UF_local->add(u, u);
		representative_local[u] = u;

		ui tu = u;
		while(true) {
			pend[tu] = pstart[tu];
			if(adj_next[tu] == tu) break;
			else tu = adj_next[tu];
		}
		adj_next[u] = u;
	}

	return c_n;
}

ui Graph::get_active_ids(const ui s, const ui *active_component, ui *active_ids, ui *degree, char *vis, ui *pend, const ui *pend_global, UnionFind *UF, const ui *representative, ui *adj_next, ui *adj_next_child) {
	// get active ids of the component of s
	// tighten adj_next, assign adj_next_child, get degree, and set pend
	// assume that all neighbors with the same component as s are put at the beginning of the adjacency list

	ui cid = active_component[s];
	ui active_n = 1;
	active_ids[0] = s; vis[s] = 1;
	for(ui i = 0;i < active_n;i ++) {
		ui u = active_ids[i];
		ui tu = u, pre = u;
		degree[u] = 0;
		while(true) {
			assert(representative[UF->UF_find(tu)] == u);
			pend[tu] = pend_global[tu];
			//if(s == 9) {
			//	printf("pstart[%u]: %u, neighbors:", tu, pstart[tu]);
			//	for(ui j = pstart[tu];j < pend_global[tu];j ++) printf(" %u", edges[j]);
			//	printf("\n");
			//}
			for(ui j = pstart[tu];j < pend_global[tu];j ++) {
				ui rv = representative[UF->UF_find(edges[j])];
				assert(rv != u);
				//if(rv == u) {
				//	swap(edges[j], edges[pstart[tu] ++]);
				//	continue;
				//}

				if(active_component[rv] != cid) {
					pend[tu] = j;
#ifndef NDEBUG
					// may not correct
					for(ui k = j;k < pend_global[tu];k ++) assert(active_component[representative[UF->UF_find(edges[k])]] != cid);
#endif
					break;
				}

				if(!vis[rv]) {
					active_ids[active_n ++] = rv;
					vis[rv] = 1;
				}
			}
			ui next = adj_next[tu];
			if(pend[tu] > pstart[tu]) {
				degree[u] += pend[tu] - pstart[tu];
				adj_next[pre] = adj_next_child[pre] = tu;
				pre = tu;
			}
			if(next == tu) break;
			else tu = next;
		}
		adj_next[pre] = adj_next_child[pre] = pre;
	}
	for(ui i = 0;i < active_n;i ++) vis[active_ids[i]] = 0;
	return active_n;
}

void Graph::read_graph_binary() {
	printf("# Start reading graph, Require files \"b_degree.bin\" and \"b_adj.bin\"\n");
	FILE *f = Utility::open_file((dir + string("/b_degree.bin")).c_str(), "rb");

	ui tt;
	fread(&tt, sizeof(ui), 1, f);
	if(tt != sizeof(ui)) {
		printf("sizeof unsigned int is different: b_degree.bin(%u), machine(%lu)\n", tt, sizeof(ui));
		return ;
	}

	fread(&n, sizeof(ui), 1, f);
	fread(&m, sizeof(ui), 1, f);

	printf("*\tn = %s; m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m/2).c_str());

	ui *degree = new ui[n];
	fread(degree, sizeof(ui), n, f);

	fclose(f);

#ifndef NDEBUG
	long long sum = 0;
	for(ui i = 0;i < n;i ++) sum += degree[i];
	assert(sum == m);
#endif

	f = Utility::open_file((dir + string("/b_adj.bin")).c_str(), "rb");

	if(pstart != nullptr) delete[] pstart;
	pstart = new ui[n+1];
	if(edges != nullptr) delete[] edges;
	edges = new ui[m];

	pstart[0] = 0;
	for(ui i = 0;i < n;i ++) {
		if(degree[i] > 0) fread(edges+pstart[i], sizeof(ui), degree[i], f);
		pstart[i+1] = pstart[i] + degree[i];
	}

	fclose(f);

	delete[] degree;

#ifndef NDEBUG
	printf("Finished reading graph\n");
#endif
}

ui Graph::core_decomposition(ui *peel_sequence, ui *core) {
#ifndef NDEBUG
	printf("Start core decomposition\n");
#endif
	ui *degree = new ui[n];
	for(ui i = 0;i < n;i ++) degree[i] = pstart[i+1]-pstart[i];

	ui *rid = new ui[n];
	ui *id = peel_sequence;
	memset(id, 0, sizeof(ui)*n);
	for(ui i = 0;i < n;i ++) ++ id[degree[i]];
	for(ui i = 1;i < n;i ++) id[i] += id[i-1];

	for(ui i = 0;i < n;i ++) rid[i] = -- id[degree[i]];
	for(ui i = 0;i < n;i ++) id[rid[i]] = i;

	ui *degree_start = new ui[n+1];
	for(ui i = 0, j = 0;i <= n;i ++) {
		while(j < n&&degree[id[j]] < i) ++ j;
		degree_start[i] = j;
	}

	ui max_core = 0;
	for(ui i = 0;i < n;i ++) {
		ui u = id[i];
		assert(degree_start[degree[u]] == i);
		if(degree[u] > max_core) max_core = degree[u];
		core[u] = max_core;

		++ degree_start[degree[u]];
		if(degree[u] == 0) continue;

		degree_start[degree[u]-1] = degree_start[degree[u]];
		for(ui j = pstart[u];j < pstart[u+1];j ++) if(rid[edges[j]] > i) {
			ui v = edges[j];
			ui pos1 = degree_start[degree[v]], pos2 = rid[v];
			swap(id[pos1], id[pos2]);
			rid[id[pos1]] = pos1; rid[id[pos2]] = pos2;
			++ degree_start[degree[v]];
			-- degree[v];
		}
	}

#ifndef NDEBUG
	for(ui i = 0;i < n;i ++) {
		ui cnt = 0;
		ui u = peel_sequence[i];
		for(ui j = pstart[u];j < pstart[u+1];j ++) if(rid[edges[j]] > i) ++ cnt;
		assert(cnt == degree[u]);
	}
#endif

	delete[] degree;
	delete[] degree_start;
	delete[] rid;

#ifndef NDEBUG
	printf("Finished core decomposition\n");
#endif
	return max_core;
}

void Graph::k_core_prune(ui K, ui *Q, ui Q_n, char *computed, ui *degree, ui *pend) {
	for(ui i = 0;i < Q_n;i ++) {
		ui u = Q[i];
		for(ui j = pstart[u];j < pend[u];j ++) if(!computed[edges[j]]) {
			ui v = edges[j];
			if(degree[v] == K) {
				Q[Q_n ++] = v;
				computed[v] = 1;
			}
			-- degree[v];
		}
	}
#ifndef NDEBUG
	printf("k_core based pruning removed %u vertices!\n", Q_n);
#endif
}

ui Graph::k_core_prune(const ui K, const ui *pend, ui *Q, ui Q_n, char *computed, ui *degree, UnionFind *UF, const ui *representative, const ui *adj_next) {
	for(ui i = 0;i < Q_n;i ++) {
		ui u = Q[i];
		ui tu = u;
		while(true) {
			assert(representative[UF->UF_find(tu)] == u);
			for(ui j = pstart[tu];j < pend[tu];j ++) {
				ui rv = representative[UF->UF_find(edges[j])];
				assert(rv != u);

				if(computed[rv]) continue;

				if(degree[rv] == K) {
					Q[Q_n ++] = rv;
					computed[rv] = 1;
				}
				-- degree[rv];
			}
			if(adj_next[tu] == tu) break;
			else tu = adj_next[tu];
		}
	}
#ifndef NDEBUG
	//printf("k_core based pruning removed %u vertices!\n", Q_n);
#endif
	return Q_n;
}

void Graph::construct_pgraph(ui s, ui *Q, char *vis, char *computed, ui *pend, Edge **graph_head, Edge *graph_edges, ui *sv_next, ui *sv_last) {
	ui cnt = 0;
	Q[0] = s; vis[s] = 1;
	ui Q_n = 1;
	assert(!computed[s]);

	for(ui i = 0;i < Q_n;i ++) {
		ui u = Q[i];
		sv_next[u] = sv_last[u] = u;
		for(ui j = pstart[u];j < pend[u];) {
			ui v = edges[j];
			if(computed[v]) swap(edges[j], edges[-- pend[u]]);
			else {
				if(!vis[v]) {
					Q[Q_n ++] = v;
					vis[v] = 1;
				}

				if(v > u) {
					graph_edges[cnt].vertex = v;
					graph_edges[cnt].reverse = &graph_edges[cnt+1];
					add_edge(u, &graph_edges[cnt], graph_head);
					++ cnt;

					graph_edges[cnt].vertex = u;
					graph_edges[cnt].reverse = &graph_edges[cnt-1];
					add_edge(v, &graph_edges[cnt], graph_head);
					++ cnt;
				}
				++ j;
			}
		}
	}
	for(ui i = 0;i < Q_n;i ++) vis[Q[i]] = 0;

	/*for(ui i = 0;i < Q_n;i ++) {
		printf("neighbors of %u:", Q[i]);
		for(Edge *e = graph_head[Q[i]];e != nullptr;e = e->next) printf(" %u", e->vertex);
		printf("\n");
	}*/
}

void Graph::initialize_pgraph(ui s,ui *Q, char *vis, char *computed, ui *pend, ui *pend_local, UnionFind *UF, ui *representative, ui *sv_next, ui *sv_last, ui *adj_next, ui *adj_last) {
	ui Q_n = 1; Q[0] = s; vis[s] = 1;
	for(ui j = 0;j < Q_n;j ++) {
		ui u = Q[j];
		for(ui k = pstart[u];k < pend[u];) {
			ui v = edges[k];
			if(computed[v]) swap(edges[k], edges[-- pend[u]]);
			else {
				if(!vis[v]) {
					Q[Q_n ++] = v;
					vis[v] = 1;
				}
				++ k;
			}
		}
	}
	UF->init(Q, Q_n);
	for(ui j = 0;j < Q_n;j ++) {
		ui u = Q[j];
		vis[u] = 0;
		representative[u] = u;
		sv_next[u] = sv_last[u] = adj_next[u] = adj_last[u] = u;
		pend_local[u] = pend[u];
	}
}

ui Graph::initialize_pgraph(const ui s, ui *pend, ui *Q, char *vis, char *computed, UnionFind *UF, const ui *representative, ui *adj_next, ui *pend_local, UnionFind *UF_local, ui *representative_local, ui *sv_next_local, ui *sv_last_local, ui *adj_next_local, ui *adj_last_local) {
	// set pend_local, adj_next_local, adj_last_local, UF_local, representative_local, sv_next_local, sv_last_local

#ifndef NDEBUG
	printf("** process a connected component\n");
#endif
	ui Q_n = 1; Q[0] = s; vis[s] = 1;
	for(ui i = 0;i < Q_n;i ++) {
		ui u = Q[i];
		ui tu = u, pre = u;
#ifndef NDEBUG
		printf("** neighbors of %u:", u);
#endif
		while(true) {
			assert(representative[UF->UF_find(tu)] == u);
			for(ui j = pstart[tu];j < pend[tu];) {
				ui rv = representative[UF->UF_find(edges[j])];
				assert(rv != u);
				if(computed[rv]) swap(edges[j], edges[-- pend[tu]]);
				else {
#ifndef NDEBUG
					printf(" %u", rv);
#endif
					UF_local->add(edges[j], rv);
					if(!vis[rv]) {
						Q[Q_n ++] = rv;
						vis[rv] = 1;
					}
					++ j;
				}
			}
			pend_local[tu] = pend[tu];
			ui next = adj_next[tu];
			if(pend[tu] > pstart[tu]) {
				adj_next[pre] = adj_next_local[pre] = tu;
				pre = tu;
			}
			if(next == tu) break;
			else tu = next;
		}
		adj_next_local[pre] = pre;
		adj_last_local[u] = pre;
#ifndef NDEBUG
		printf("\n");
#endif
	}
	UF_local->init(Q, Q_n);

	for(ui j = 0;j < Q_n;j ++) {
		ui u = Q[j];
		vis[u] = 0;
		representative_local[u] = u;
		sv_next_local[u] = sv_last_local[u] = u;
	}

	return Q_n;
}

ui Graph::decomposition(ui s, ui K, ui *cstart, ui *ids, ui c_n, ui *Q, ui *keys, char *vis, Edge **graph_head, Edge *graph_edges, ListLinearHeap *heap, ui *sv_next, ui *sv_last) {
	while(graph_head[s] != nullptr) {
		heap->init(0, K, nullptr, nullptr);
		heap->insert(s, 0);
		ui u, key, Q_n = 0;
		while(heap->pop_max(u, key)) {
			Q[Q_n ++] = u; vis[u] = 1; //vis[u] = 1 means u is in Q, and vis[u] = 2 means u is in heap
										//vis[u] = 3 means u is in Q between Q_n and new_Qn
			keys[u] = key;

			ui new_Qn = Q_n;
			for(ui j = Q_n-1;j < new_Qn;j ++) {
				ui v = Q[j];

				for(Edge *e = graph_head[v];e != nullptr;e = e->next) if(vis[e->vertex] != 1) {
					ui w = e->vertex;

					if(vis[w] == 3) {
						++ keys[w];
						continue;
					}

					if(vis[w] == 2) key = heap->remove(w);
					else key = 0;
					assert(key < K);

					++ key;
					if(key >= K) {
						Q[new_Qn ++] = w;
						keys[w] = key;
						vis[w] = 3;
					}
					else {
						heap->insert(w, key);
						vis[w] = 2;
					}
				}

				if(v == u) continue;

				// contract u and v
				vis[v] = 0;
				keys[u] += keys[v];
				merge(graph_head, graph_edges, u, v, keys, sv_next, sv_last);
			}
		}

		/*printf("before\n");
		for(ui i = 0;i < Q_n;i ++) {
			printf("neigbors of %u:", Q[i]);
			for(Edge *e = graph_head[Q[i]];e != nullptr;e = e->next) printf(" %u", e->vertex);
			printf("\n");
		}*/

		while(Q_n > 0&&keys[Q[Q_n-1]] < K) {
			ui u = Q[-- Q_n];
			vis[u] = 0;
			//printf("after remove %u\n", u);

			for(Edge *e = graph_head[u];e != nullptr;e = e->next) delete_edge(e->vertex, e->reverse, graph_head);
			graph_head[u] = nullptr;

			ui pos = cstart[c_n];
			ids[pos ++] = u;
			while(sv_next[u] != u) {
				u = sv_next[u];
				ids[pos ++] = u;
			}
			cstart[++ c_n] = pos;

			/*for(ui i = 0;i < Q_n;i ++) {
				printf("neigbors of %u:", Q[i]);
				for(Edge *e = graph_head[Q[i]];e != nullptr;e = e->next) printf(" %u", e->vertex);
				printf("\n");
			}*/
		}

		for(ui i = 0;i < Q_n;i ++) vis[Q[i]] = 0;
	}
	return c_n;
}

ui Graph::decomposition(const ui s, const ui K, ui *cstart, ui *ids, ui c_n, ui *Q, ui *keys, char *vis, ui *pend_local, ListLinearHeap *heap, ui *sv_next, ui *sv_last, ui *adj_next, ui *adj_last, UnionFind *UF, ui *representative) {
	ui old_cn = c_n;
	while(true) {
		heap->init(0, K, nullptr, nullptr);
		heap->insert(s, 0);
		ui u, key, Q_n = 0;
		//printf("start an iteration\n");
		while(heap->pop_max(u, key)) {
			Q[Q_n ++] = u; vis[u] = 1; //vis[u] = 1 means u is in Q, and vis[u] = 2 means u is in heap
										//vis[u] = 3 means u is in Q between Q_n and new_Qn
										//vis[u] = 4 means u is removed from the pgraph
			keys[u] = key;
			//printf(" [%u,%u]\n", u, key);

			ui new_Qn = Q_n;
			for(ui i = Q_n-1;i < new_Qn;i ++) {
				ui v = Q[i];
				ui pre = v, tv = v;

				while(true) {
					assert(representative[UF->UF_find(v)] == v);
					for(ui j = pstart[tv];j < pend_local[tv];) {
						ui w = representative[UF->UF_find(edges[j])];
						if(vis[w] == 4||w == v) {
							swap(edges[j], edges[-- pend_local[tv]]);
							continue;
						}
						++ j;
						if(vis[w] != 1) {
							if(vis[w] == 3) {
								++ keys[w];
								continue;
							}

							if(vis[w] == 2) key = heap->remove(w);
							else key = 0;
							assert(key < K);

							++ key;
							if(key >= K) {
								Q[new_Qn ++] = w;
								keys[w] = key;
								vis[w] = 3;
							}
							else {
								heap->insert(w, key);
								vis[w] = 2;
							}
						}
					}
					ui next = adj_next[tv];
					if(pend_local[tv] > pstart[tv]) {
						adj_next[pre] = tv;
						pre = tv;
					}
					if(next == tv) break;
					else tv = next;
				}
				adj_last[v] = pre;
				adj_next[pre] = pre;

				if(v == u) continue;

				// contract u and v
				keys[u] += keys[v];
				//printf("merged %u and %u\n", u, v);
				tv = v;
				while(true) {
					for(ui j = pstart[tv];j < pend_local[tv];j ++) if(representative[UF->UF_find(edges[j])] == u) -- keys[u];
					if(adj_next[tv] == tv) break;
					else tv = adj_next[tv];
				}
				sv_next[sv_last[u]] = v;
				sv_last[u] = sv_last[v];
				adj_next[adj_last[u]] = v;
				adj_last[u] = adj_last[v];
				representative[UF->UF_union(u,v)] = u;
			}
		}

		while(Q_n > 0&&keys[Q[Q_n-1]] < K) {
			ui u = Q[-- Q_n];

			ui pos = cstart[c_n];
			ids[pos ++] = u; vis[u] = 4;
			while(sv_next[u] != u) {
				u = sv_next[u];
				ids[pos ++] = u; vis[u] = 4;
			}
			cstart[++ c_n] = pos;
		}

		for(ui i = 0;i < Q_n;i ++) vis[Q[i]] = 0;
		if(Q_n == 0) break;
		//printf("finished one iteration\n\n");
	}
	for(ui i = old_cn;i < c_n;i ++) for(ui j = cstart[i];j < cstart[i+1];j ++) vis[ids[j]] = 0;

	return c_n;
}

void Graph::merge(Edge **graph_head, Edge *edges, ui u, ui v, ui *keys, ui *sv_next, ui *sv_last) {
	for(Edge *e = graph_head[v];e != nullptr;) {
		Edge *tmp = e->next;

		if(e->vertex == u) {
			-- keys[u];
			delete_edge(u, e->reverse, graph_head);
		}
		else {
			assert(e->reverse->vertex == v);
			e->reverse->vertex = u;

			add_edge(u, e, graph_head);
		}

		e = tmp;
	}
	graph_head[v] = nullptr;

	sv_next[sv_last[u]] = v;
	sv_last[u] = sv_last[v];
}

void Graph::delete_edge(ui u, Edge *e, Edge **graph_head) {
	if(e->pre == nullptr) {
		assert(graph_head[u] == e);
		e = e->next;
		if(e != nullptr) e->pre = nullptr;
		graph_head[u] = e;
	}
	else {
		assert(graph_head[u] != e);
		e->pre->next = e->next;
		if(e->next != nullptr) e->next->pre = e->pre;
	}
}

void Graph::add_edge(ui u, Edge *e, Edge **graph_head) {
	if(graph_head[u] != nullptr) graph_head[u]->pre = e;
	e->next = graph_head[u];
	graph_head[u] = e;
	e->pre = nullptr;
}

void Graph::remove_inter_edges(ui K, ui c_n, ui new_cn, ui *cstart, ui *ids, ui *component_id, ui *pend, ui *Q, char *computed, ui *degree) {
	for(ui i = c_n;i < new_cn;i ++) for(ui j = cstart[i];j < cstart[i+1];j ++) component_id[ids[j]] = i;

	for(ui i = c_n;i < new_cn;i ++) for(ui j = cstart[i];j < cstart[i+1];j ++) {
		ui u = ids[j];
		for(ui k = pstart[u];k < pend[u];) {
			ui v = edges[k];
			assert(!computed[v]);
			if(component_id[v] != component_id[u]) swap(edges[k], edges[-- pend[u]]);
			else ++ k;
		}
	}

	ui Q_n = 0;
	for(ui i = c_n;i < new_cn;i ++) for(ui j = cstart[i];j < cstart[i+1];j ++) {
		ui u = ids[j];
		degree[u] = pend[u] - pstart[u];
		if(degree[u] < K) {
			Q[Q_n ++] = u;
			computed[u] = 1;
		}
	}

	k_core_prune(K, Q, Q_n, computed, degree, pend) ;
}

void Graph::remove_inter_edges(ui K, ui c_n, ui new_cn, ui *cstart, ui *ids, ui *pend, ui *Q, char *computed, ui *degree, UnionFind *UF) {
	ui Q_n = 0;
	for(ui i = c_n;i < new_cn;i ++) for(ui j = cstart[i];j < cstart[i+1];j ++) {
		ui u = ids[j];
		ui ru = UF->UF_find(u);
		for(ui k = pstart[u];k < pend[u];) {
			if(UF->UF_find(edges[k]) != ru) swap(edges[k], edges[-- pend[u]]);
			else ++ k;
		}
		degree[u] = pend[u] - pstart[u];
		if(degree[u] < K) {
			Q[Q_n ++] = u;
			computed[u] = 1;
		}
	}
	k_core_prune(K, Q, Q_n, computed, degree, pend);
}

ui Graph::remove_inter_edges(const ui K, const ui c_n, const ui new_cn, const ui *cstart, const ui *ids, ui *pend, ui *Q, char *computed, ui *degree, UnionFind *UF_local, const ui *representative_local, ui *adj_next) {
	ui Q_n = 0;
	for(ui i = c_n;i < new_cn;i ++) for(ui j = cstart[i];j < cstart[i+1];j ++) {
		ui u = ids[j];
		ui tu = u, pre = u, ru = UF_local->UF_find(u);
		degree[u] = 0;
		while(true) {
			for(ui k = pstart[tu];k < pend[tu];) {
				if(UF_local->UF_find(edges[k]) != ru) swap(edges[k], edges[-- pend[tu]]);
				else ++ k;
			}

			ui next = adj_next[tu];
			if(pend[tu] > pstart[tu]) {
				degree[u] += pend[tu] - pstart[tu];
				adj_next[pre] = tu;
				pre = tu;
			}
			if(next == tu) break;
			else tu = next;
		}
		adj_next[pre] = pre;

		if(degree[u] < K) {
			assert(!computed[u]);
			Q[Q_n ++] = u;
			computed[u] = 1;
		}
	}
	return Q_n;
}

void Graph::print_keccs(ui level, ui K, ui c_n, ui *cstart, ui *ids) {
	for(ui i = 0;i < c_n;i ++) {
		for(ui j = 0;j < level;j ++) printf("*");
		printf(" %u-ecc %u:", K, i+1);
		for(ui k = cstart[i];k < cstart[i+1];k ++) printf(" %u", ids[k]);
		printf("\n");
	}
}

void Graph::print_active_graph(const ui *active_ids, const ui active_n, const ui *pend, const ui *adj_next, UnionFind *UF, const ui *representative) {
	printf("* Active ids:");
	for(ui i = 0;i < active_n;i ++) printf(" %u", active_ids[i]);
	printf("\n");

	for(ui i = 0;i < active_n;i ++) {
		ui u = active_ids[i];
		ui tu = u;
		printf("* Neigbors of %u:", u);
		while(true) {
			assert(representative[UF->UF_find(tu)] == u);
			for(ui j = pstart[tu];j < pend[tu];j ++) {
				ui rv = representative[UF->UF_find(edges[j])];
				assert(rv != u);
				printf(" %u", rv);
			}
			if(adj_next[tu] == tu) break;
			else tu = adj_next[tu];
		}
		printf("\n");
	}
	printf("\n");
}

void Graph::print_kecc(string output_file, ui K, ui c_n, ui *cstart, ui *ids) {
	//ostringstream os;
	//os<<K<<"-ECCs.txt";
	printf("Writing kECCS to file: %s\n", output_file.c_str());
	FILE *fout = Utility::open_file(output_file.c_str(), "w");

	for(ui i = 0;i < c_n;i ++) {
		fprintf(fout, "%u-ECC %u of size %u:", K, i, cstart[i+1]-cstart[i]);
		sort(ids+cstart[i], ids+cstart[i+1]);
		for(ui j = cstart[i];j < cstart[i+1];j ++) fprintf(fout, " %u", ids[j]);
		fprintf(fout, "\n");
	}

	fclose(fout);
}

void Graph::print_edge_connectivities(const char *alg, vector<pair<pair<ui,ui>, ui> > &vpp, bool mspt) {
	ostringstream os;
	if(mspt) os<<alg<<"-mspt.txt";
	else os<<alg<<"-all.txt";
	FILE *fout = Utility::open_file(os.str().c_str(), "w");

	fprintf(fout, "u v sc\n");
	sort(vpp.begin(), vpp.end());
	printf("*\tvpp has been sorted based on vertex ids!!!\n");
	for(ui i = 0;i < vpp.size();i ++) fprintf(fout, "%u %u %u\n", vpp[i].first.first, vpp[i].first.second, vpp[i].second);

	fclose(fout);
}
