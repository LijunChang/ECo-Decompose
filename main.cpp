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
	//print_usage();

	OptionParser op("Allowed options");
	auto help_option = op.add<Switch>("h", "help", "\'produce help message\'");
	auto graph_option = op.add<Value<string>>("g", "graph", "\'path to input graph file\'");
	auto alg_option = op.add<Value<string>>("a", "alg", "\'algorithm name\' (kecc-space | eco-decompose-dcs)");
	auto k_option = op.add<Value<int>>("k", "k", "\'the value of k for kecc-space\'");
	auto output_option = op.add<Value<string>>("o", "output", "\'write the result into file\'");
	auto eccsizes_option = op.add<Value<string>>("s", "eccsizes", "\'write the ecc_sizes into file\'");

	op.parse(argc, argv);

	if(help_option->is_set()||argc <= 1) {
		cout << op << endl;
		if(argc <= 1) {
			print_usage();
			return 0;
		}
	}
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

	if(strcmp(alg.c_str(), "kecc-space") == 0) {
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

	string output_file = "", eccsizes_file = "";
	if(output_option->is_set()) output_file = output_option->value();
	if(eccsizes_option->is_set()) eccsizes_file = eccsizes_option->value();

	Graph *graph = new Graph(graph_option->value().c_str());
	graph->read_graph_binary();

	Timer timer;
	if(strcmp(alg.c_str(), "kecc-space") == 0) graph->k_edge_connected_component_space((ui)k, output_file);
	else if(strcmp(alg.c_str(), "eco-decompose-dcs") == 0) graph->edge_connectivity_decomposition_DCs(true, output_file, eccsizes_file);
	else {
		printf("!!! The algorithm name is not reconganized! Exit!!!\n");
		//print_usage();
	}
	delete graph;

	printf("* Total processing time excluding I/O: %s (microseconds)\n", Utility::integer_to_string(timer.elapsed()).c_str());
	if(!output_file.empty()||!eccsizes_file.empty()) printf("*\tNote that this includes the time of writing the output to file!\n");
	printf("**********************************\n");

	return 0;
}

void print_usage() {
	printf("Example usage: ./eco_decompose -g path_to_graph -a kecc-space -k 3\n");
	printf("or\t./eco_decompose -g path_to_graph -a eco-decompose-dcs\n");
}
