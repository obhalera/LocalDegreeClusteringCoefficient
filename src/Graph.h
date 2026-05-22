#include <iostream>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <algorithm>
#include <random>
#include <stdexcept>
#include <chrono>  
using namespace std;

class CSRGraph {
private:
    long long num_nodes = 0;
    long long num_edges = 0;

    vector<long long> row_ptr;
    vector<long long> col_idx;

 
    mt19937 rng;
public:
    
    // ============================================================
    // CONSTRUCTOR: robust CSR builder (directed/undirected safe)
    // ============================================================
    CSRGraph(const string& filepath, bool directed_input = false) {

        ifstream infile(filepath);
        if (!infile)
            throw runtime_error("Could not open file");

        vector<pair<long long,long long>> edges;
        edges.reserve(1e6);

        unordered_map<long long,long long> id_map;
        id_map.reserve(1e6);

        long long u, v;
        long long next_id = 0;

        // -----------------------------
        // Step 1: read + compress IDs
        // -----------------------------
        while (infile >> u >> v) {

            if (!id_map.count(u)) id_map[u] = next_id++;
            if (!id_map.count(v)) id_map[v] = next_id++;

            long long a = id_map[u];
            long long b = id_map[v];

            if (a == b) continue; // remove self-loops

            edges.emplace_back(a, b);

            // symmetrize if directed input
            if (directed_input)
                edges.emplace_back(b, a);
        }

        num_nodes = next_id;

        // -----------------------------
        // Step 2: deduplicate edges
        // -----------------------------
        sort(edges.begin(), edges.end());
        edges.erase(unique(edges.begin(), edges.end()), edges.end());

        num_edges = edges.size();

        // -----------------------------
        // Step 3: build CSR structure
        // -----------------------------
        row_ptr.assign(num_nodes + 1, 0);

        for (auto &e : edges)
            row_ptr[e.first + 1]++;

        for (long long i = 1; i <= num_nodes; i++)
            row_ptr[i] += row_ptr[i - 1];

        col_idx.assign(num_edges, 0);

        vector<long long> temp_ptr = row_ptr;

        for (auto &e : edges)
            col_idx[temp_ptr[e.first]++] = e.second;

        // -----------------------------
        // Step 4: sort adjacency lists
        // -----------------------------
        for (long long u = 0; u < num_nodes; u++) {
            sort(col_idx.begin() + row_ptr[u],
                 col_idx.begin() + row_ptr[u + 1]);
        }

        // -----------------------------
        // RNG init
        // -----------------------------
        random_device rd;
        rng.seed(rd() ^ chrono::high_resolution_clock::now()
                    .time_since_epoch().count());
    }

    // ============================================================
    // BASIC GRAPH PROPERTIES
    // ============================================================
    long long numVertices() const {
        return num_nodes;
    }

    long long numEdges() const {
        return num_edges;
    }

    long long degree(long long u) const {
        if (u < 0 || u >= num_nodes)
            throw out_of_range("Invalid vertex");

        return row_ptr[u + 1] - row_ptr[u];
    }

    // ============================================================
    // NEIGHBOR ACCESS (FAST + DEBUG)
    // ============================================================
    vector<long long> neighbors(long long u) const {
        return vector<long long>(
            col_idx.begin() + row_ptr[u],
            col_idx.begin() + row_ptr[u + 1]
        );
    }

    inline const long long* neighbors_begin(long long u) const {
        return &col_idx[row_ptr[u]];
    }

    inline const long long* neighbors_end(long long u) const {
        return &col_idx[row_ptr[u + 1]];
    }

    // ============================================================
    // EDGE QUERY (BINARY SEARCH)
    // ============================================================
    bool has_edge(long long u, long long v) const {
        if (u < 0 || u >= num_nodes) return false;

        return binary_search(
            col_idx.begin() + row_ptr[u],
            col_idx.begin() + row_ptr[u + 1],
            v
        );
    }

    // ============================================================
    // RANDOM NEIGHBOR (UNIFORM)
    // ============================================================
    long long random_neighbor(long long u) {
        long long deg = degree(u);
        if (deg == 0)
            throw runtime_error("Isolated vertex");

        uniform_int_distribution<long long> dist(0, deg - 1);
        return col_idx[row_ptr[u] + dist(rng)];
    }

    std::vector<double> actual_degree_clustering_coefficient();
    std::vector<double> estimated_degree_clustering_coefficient(double frac_nodes, double frac_edges);
    std::vector<double> edge_estimated_degree_clustering_coefficient(double r_frac_edges, double l_frac_edges);
    void compare_clustering_coefficients(const std::vector<double>& exact, const std::vector<double>& estimated, const std::string& output_file);
 
    std::vector<long long> random_walk(long long start, long long steps);

    std::vector<double> rw_vertex_estimated_degree_clustering_coefficient(double r_frac_edges, double l_frac_edges, vector<long long> &walk_vertices);
    std::vector<double> clustering_coeff_sparsification(double p, vector<long long> &walk);

    std::vector<double> rw_degree_cycle_estimation(double r_frac_edges, double l_frac_edges, vector<long long> &red_nodes, vector<long long> &walk_vertices);
    std::vector<double> degree_4cycle_distribution(vector<long long> &red_vertices);
};