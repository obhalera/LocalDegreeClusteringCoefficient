#include <iostream>
#include <vector>
#include <string>
#include "Graph.h"
//Compile: clang++ main.cpp graph.cpp -std=c++17
//./main graph.txt 0.1 0.1
int main(int argc, char* argv[]) {

    // Read graph file path from command line argument
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <graph_file(without ext)> [frac_nodes] [frac_edges] \n";
        return 1;
    }

    std::string input_base_path = "/Users/omkarbhalerao/Desktop/LocalDegreeClusteringCoefficient/datasets/";
    std::string output_base_path = "/Users/omkarbhalerao/Desktop/LocalDegreeClusteringCoefficient/output/distributions/";
    std::string curr_graph_string;
    std::string dataset = argv[1];
    std::string filepath = input_base_path + argv[1] + ".txt";
    std::string curr_method = "estimators.csv";


    if(dataset == "facebook")
        curr_graph_string = "facebook/";
    else if(dataset == "com-youtube-sanitized")
        curr_graph_string = "youtube/";
    else if(dataset == "soc-LiveJournal1-sanitized")
        curr_graph_string = "livejournal/";
    else if(dataset == "soc-orkut-sanitized")
        curr_graph_string = "orkut/";
    else if(dataset == "wiki-topcats-sanitized")
        curr_graph_string = "wiki-topcats/";

    std::string outputpath = output_base_path + curr_graph_string + curr_method;

    // Optional: read frac_nodes and frac_edges from command line, default to 0.1
    double frac_nodes = (argc >= 2) ? std::stod(argv[2]) : 0.1;
    double frac_edges = (argc >= 3) ? std::stod(argv[3]) : 0.1;

    

    // Load graph
    std::cout << "Loading graph from: " << filepath << "\n";
    CSRGraph G(filepath, /*undirected=*/true);
    std::cout << "Nodes: " << G.numVertices()
              << "  Edges: " << G.numEdges() << "\n\n";

    long long steps = (long long)((frac_nodes + frac_edges) * (G.numEdges()));

    // Step 1: Exact degree clustering coefficients
    std::cout << "Computing exact degree clustering coefficients...\n";
    std::vector<double> exact = G.actual_degree_clustering_coefficient();

    //for(auto i = 0; i < exact.size(); i++)
    //    std::cout<<"Degree: "<<i<<", Returned triangle: "<<exact[i]<<endl;
    std::cout << "Done.\n\n";

    mt19937 rng_seed_sampler;
    random_device rd_seed_sampler;
    rng_seed_sampler.seed(rd_seed_sampler() ^ chrono::high_resolution_clock::now()
                    .time_since_epoch().count());

    std::uniform_int_distribution<long long> node_dist(0, G.numVertices() - 1);
    long long seed = node_dist(rng_seed_sampler);

    vector<long long> walk_vertices = G.random_walk(seed, steps);
    vector<long long> vtx_est_walk_vertices;

    for(auto i = 0; i < (long long)(frac_nodes * walk_vertices.size()); i++) vtx_est_walk_vertices.push_back(walk_vertices[i]);

    std::cout<<"Length of the original random walk: "<<walk_vertices.size()<<endl;
    std::cout<<"Length of the walk fed as input to rw walk estimator: "<<vtx_est_walk_vertices.size()<<endl;

    std::cout<<"Computing the estimated clustering coefficient using edge samples...\n";
    std::vector<double> estimated = G.rw_vertex_estimated_degree_clustering_coefficient(frac_nodes, frac_edges, vtx_est_walk_vertices);
    std::cout<<"Done.\n\n";

    // Step 3: Compare and print relative errors per bin
    std::cout << "Comparing exact vs estimated:\n\n";
    G.compare_clustering_coefficients(exact, estimated, outputpath);    //Later need to update this function for different baselines

    return 0;
}