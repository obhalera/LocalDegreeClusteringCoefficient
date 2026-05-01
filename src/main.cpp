#include <iostream>
#include <vector>
#include <string>
#include "Graph.h"
//Compile: clang++ main.cpp graph.cpp -std=c++17
//./main graph.txt 0.1 0.1
int main(int argc, char* argv[]) {

    // Read graph file path from command line argument
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <graph_file> [frac_nodes] [frac_edges]\n";
        return 1;
    }

    std::string filepath = argv[1];

    // Optional: read frac_nodes and frac_edges from command line, default to 0.1
    double frac_nodes = (argc >= 3) ? std::stod(argv[2]) : 0.1;
    double frac_edges = (argc >= 4) ? std::stod(argv[3]) : 0.1;

    // Load graph
    std::cout << "Loading graph from: " << filepath << "\n";
    CSRGraph G(filepath, /*undirected=*/true);
    std::cout << "Nodes: " << G.numVertices()
              << "  Edges: " << G.numEdges() << "\n\n";

    // Step 1: Exact degree clustering coefficients
    std::cout << "Computing exact degree clustering coefficients...\n";
    std::vector<double> exact = G.actual_degree_clustering_coefficient();
    std::cout << "Done.\n\n";

    // Step 2: Estimated degree clustering coefficients
   /* std::cout << "Computing estimated degree clustering coefficients "
              << "(frac_nodes=" << frac_nodes
              << ", frac_edges=" << frac_edges << ")...\n";
    std::vector<double> estimated = G.estimated_degree_clustering_coefficient(frac_nodes, frac_edges);
    std::cout << "Done.\n\n";*/

    std::cout<<"Computing the estimated clustering coefficient using edge samples...\n";
    std::vector<double> estimated = G.edge_estimated_degree_clustering_coefficient(frac_nodes, frac_edges);
    std::cout<<"Done.\n\n";

    // Step 3: Compare and print relative errors per bin
    std::cout << "Comparing exact vs estimated:\n\n";
    G.compare_clustering_coefficients(exact, estimated);

    return 0;
}