#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <random>
#include <stdexcept>

class CSRGraph {
private:
    long long num_nodes;
    long long num_edges;

    std::vector<long long> row_ptr;   // size = num_nodes + 1
    std::vector<long long> col_idx;   // size = num_edges

    std::mt19937 rng;

public:
    // Constructor
    CSRGraph(const std::string& filepath, bool undirected = true) {

        std::random_device rd;
        auto seed = rd()
                    ^ std::chrono::high_resolution_clock::now().time_since_epoch().count();
        std::mt19937 gen(seed);

        std::ifstream infile(filepath);
        if (!infile) {
            throw std::runtime_error("Could not open file");
        }

        std::vector<std::pair<long long, long long> > edges;
        long long u, v;
        long long max_node = -1;

        // Step 1: Read edges + symmetrize if undirected
        while (infile >> u >> v) {
            //edges.emplace_back(u, v);
            edges.push_back(std::make_pair(u, v));
            if(u < v) max_node = v;
            else max_node = u;
            //max_node = std::max({max_node, u, v});

            if (undirected && u != v) {
                //edges.emplace_back(v, u);
                edges.push_back(std::make_pair(v, u));
            }
        }

        // Step 2: Deduplicate edges
        std::sort(edges.begin(), edges.end());
        edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

        num_nodes = max_node + 1;
        num_edges = edges.size();

        // Step 3: Build row_ptr (degree counts)
        row_ptr.assign(num_nodes + 1, 0);
        for (const auto& e : edges) {
            row_ptr[e.first + 1]++;
        }

        // Step 4: Prefix sum
        for (long long i = 1; i <= num_nodes; i++) {
            row_ptr[i] += row_ptr[i - 1];
        }

        // Step 5: Fill col_idx
        col_idx.assign(num_edges, 0);
        std::vector<long long> temp_ptr = row_ptr;

        for (const auto& e : edges) {
            long long src = e.first;
            long long dst = e.second;
            col_idx[temp_ptr[src]++] = dst;
        }

        // Step 6: Sort adjacency lists (for binary search)
        for (long long i = 0; i < num_nodes; i++) {
            std::sort(col_idx.begin() + row_ptr[i],
                      col_idx.begin() + row_ptr[i + 1]);
        }
    }

    // Number of vertices
    long long numVertices() const {
        return num_nodes;
    }

    // Number of edges (directed count; undirected = 2 * unique edges)
    long long numEdges() const {
        return num_edges;
    }

    // Degree of vertex u
    long long degree(long long u) const {
        if (u < 0 || u >= num_nodes) {
            throw std::out_of_range("Invalid vertex");
        }
        return row_ptr[u + 1] - row_ptr[u];
    }

    // Check if edge (u, v) exists
    bool has_edge(long long u, long long v) const {
        if (u < 0 || u >= num_nodes) return false;

        long long start = row_ptr[u];
        long long end = row_ptr[u + 1];

        return std::binary_search(
            col_idx.begin() + start,
            col_idx.begin() + end,
            v
        );
    }

    // Sample a uniform random neighbor of u
    long long random_neighbor(long long u) {
        long long deg = degree(u);
        if (deg == 0) {
            throw std::runtime_error("Vertex has no neighbors");
        }

        std::uniform_int_distribution<long long> dist(0, deg - 1);
        long long idx = dist(rng);

        return col_idx[row_ptr[u] + idx];
    }

    // Get all neighbors of u
    std::vector<long long> neighbors(long long u) const {
        if (u < 0 || u >= num_nodes) {
            throw std::out_of_range("Invalid vertex");
        }

        return std::vector<long long>(
            col_idx.begin() + row_ptr[u],
            col_idx.begin() + row_ptr[u + 1]
        );
    }

    std::vector<double> actual_degree_clustering_coefficient();
    std::vector<double> estimated_degree_clustering_coefficient(double frac_nodes, double frac_edges);
    std::vector<double> edge_estimated_degree_clustering_coefficient(double r_frac_edges, double l_frac_edges);
    void compare_clustering_coefficients(const std::vector<double>& exact, const std::vector<double>& estimated);
};