/*
 * main.cpp
 *
 *  Created on: 5Dec.,2017
 *      Author: Lijun Chang
 *      Email: ljchang@outlook.com
 */

#include "Graph.h"
#include "popl.hpp"

using namespace std;
using namespace popl;

void print_usage() ;

int main(int argc, char *argv[]) {
#ifndef NDEBUG
	printf("**** eco_decompose (Debug) build at %s %s ***\n", __TIME__, __DATE__);
#else
	printf("**** eco_decompose (Release) build at %s %s ***\n", __TIME__, __DATE__);
#endif

#ifndef NDEBUG
	printf("!!! You may want to define NDEBUG in utilities/Defines.h to get better performance!\n");
#endif
	//printf("sizeof(unsigned int): %lu\n", sizeof(ui));
	print_usage();

	OptionParser op("Allowed options");
	auto help_option = op.add<Switch>("h", "help", "\'produce help message\'");
	auto graph_option = op.add<Value<string>>("g", "graph", "\'path to input graph file\'");
	auto alg_option = op.add<Value<string>>("a", "alg", "\'algorithm name\' (kecc | kecc-space | eco-decompose-dcs | eco-decompose-bus | eco-decompose-tds)");
	auto k_option = op.add<Value<int>>("k", "k", "\'the value of k\'");
	auto output_option = op.add<Switch>("o", "output", "\'write the result into a file\'");

	op.parse(argc, argv);

	if(help_option->is_set()) cout << op << endl;
	if(!graph_option->is_set()) {
		printf("!!! Path to input graph file is not provided! Exit !!!\n");
		return 0;
	}
	if(!alg_option->is_set()) {
		printf("!!! Algorithm name is not provided! Exit !!!\n");
		return 0;
	}

	string alg = alg_option->value();
	int k = 0;

	if(strcmp(alg.c_str(), "kecc")==0 || strcmp(alg.c_str(), "kecc-space") == 0) {
		if(!k_option->is_set()) {
			printf("!!! The value of k is not provided for KECC computation! Exit !!!\n");
			return 0;
		}
		else {
			k = k_option->value();
			if(k < 2) {
				printf("!!! k must be at least 2! Exit !!!\n");
				return 0;
			}
		}
	}

	bool output = false;
	if(output_option->is_set()) output = true;

	Graph *graph = new Graph(graph_option->value().c_str());
	graph->read_graph_binary();

	Timer timer;
	if(strcmp(alg.c_str(), "kecc") == 0) graph->k_edge_connected_component((ui)k, output);
	else if(strcmp(alg.c_str(), "kecc-space") == 0) graph->k_edge_connected_component_space((ui)k, output);
	else if(strcmp(alg.c_str(), "eco-decompose-dcs") == 0) {
	//	bool mspt = false;
	//	if(argc >= 4&&strcmp(argv[3], "mspt") == 0) mspt = true;
		graph->edge_connectivity_decomposition_DCs(false, output);
	}
	//else if(strcmp(alg.c_str(), "eco-decompose-buso") == 0) {
	//	bool mspt = false;
	//	if(argc >= 4&&strcmp(argv[3], "mspt") == 0) mspt = true;
	//	graph->edge_connectivity_decomposition_BUso(mspt, print);
	//}
	else if(strcmp(alg.c_str(), "eco-decompose-bus") == 0) {
	//	bool mspt = false;
	//	if(argc >= 4&&strcmp(argv[3], "mspt") == 0) mspt = true;
		graph->edge_connectivity_decomposition_BUs(true, output);
	}
	else if(strcmp(alg.c_str(), "eco-decompose-tds") == 0) {
	//	bool mspt = false;
	//	if(argc >= 4&&strcmp(argv[3], "mspt") == 0) mspt = true;
		graph->edge_connectivity_decomposition_TDs(true, output);
	}
	else {
		printf("!!! The algorithm name is not reconganized! Exit!!!\n");
		//print_usage();
	}
	delete graph;

	printf("*\tTotal processing time excluding I/O: %s (microseconds)\n", Utility::integer_to_string(timer.elapsed()).c_str());
	if(output) printf("*\tNote that this includes the time of writing the output to file!\n");
	printf("**********************************\n");

	return 0;
}

void print_usage() {
	printf("Usage: ./eco_decompose -h -g path_to_graph -a kecc-space -k 3\n");
	printf("or\t./eco_decompose -h -g path_to_graph -a eco-decompose-dcs\n");
}
