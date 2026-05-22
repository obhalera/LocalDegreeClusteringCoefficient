#include <iostream>
#include <vector>
#include <fstream>
#include <set>
#include <map>
#include <unordered_set>
#include <algorithm>
#include <random>
#include <stdexcept>
#include "Graph.h"


//Need to update this function
std::vector<double> CSRGraph::actual_degree_clustering_coefficient() {
    
    std::vector<long long> triangle_count(num_nodes, 0);
    std::vector<long long> node_degree(num_nodes, 0);

    //for(auto i = 0; i < num_nodes; i++)
    //    std::cout<<"Node: "<<i<<", degree: "<<degree(i)<<std::endl;

    // Step 1: Compute degrees
    for (long long i = 0; i < num_nodes; i++) {
        node_degree[i] = degree(i);
    }

    // Step 2: Define total ordering (degree, then id)
    auto less_than = [&](long long a, long long b) {
        if (node_degree[a] != node_degree[b])
            return node_degree[a] < node_degree[b];
        return a < b;
    };

    // Step 3: Build forward adjacency lists
    std::vector<std::vector<long long>> forward(num_nodes);
    for (long long u = 0; u < num_nodes; u++) {
        for (long long v : neighbors(u)) {
            if (less_than(u, v)) {
                forward[u].push_back(v);
            }
        }
        std::sort(forward[u].begin(), forward[u].end());
    }

    // Step 4: Triangle counting via intersection
    for (long long u = 0; u < num_nodes; u++) {
        for (long long v : forward[u]) {
            long long i = 0, j = 0;

            const auto& Nu = forward[u];
            const auto& Nv = forward[v];

            while (i < (long long)Nu.size() && j < (long long)Nv.size()) {
                if (Nu[i] == Nv[j]) {
                    long long w = Nu[i];

                    // triangle (u, v, w) counted exactly once
                    triangle_count[u]++;
                    triangle_count[v]++;
                    triangle_count[w]++;

                    i++;
                    j++;
                } else if (Nu[i] < Nv[j]) {
                    i++;
                } else {
                    j++;
                }
            }
        }
    }

    // Step 5: Max degree
    long long max_degree = *std::max_element(node_degree.begin(), node_degree.end());

    // Step 6: Count nodes per degree
    std::vector<long long> degree_node_count(max_degree + 1, 0);
    for (long long i = 0; i < num_nodes; i++) {
        degree_node_count[node_degree[i]]++;
    }

    // Step 7: Sum triangle counts per degree
    std::vector<double> degree_triangle_sum(max_degree + 1, 0.0);
    for (long long i = 0; i < num_nodes; i++) {
        degree_triangle_sum[node_degree[i]] += (double)triangle_count[i];
    }

    //for(auto i = 0; i < degree_triangle_sum.size(); i++)
    //    std::cout<<"Degree: "<<i<<", Num triangles: "<<degree_triangle_sum[i]<<std::endl;

    // Step 8: Compute clustering (you had this partially disabled)
    std::vector<double> clustering(max_degree + 1, 0.0);
    for (long long d = 2; d <= max_degree; d++) {
        if (degree_node_count[d] == 0) continue;

        double d_choose_2 = (double)d * (d - 1) / 2.0;

        // If you want ACTUAL clustering coefficient, use this:
        clustering[d] = degree_triangle_sum[d] /
                        (degree_node_count[d] * d_choose_2);

        // If you only want raw triangle sums, use:
        //clustering[d] = degree_triangle_sum[d];
    }

    //for(auto i = 0; i < clustering.size(); i++)
    //    cout<<"Degree: "<<i<<", clustering: "<<clustering[i]<<endl;

    return clustering;
}

std::vector<long long> CSRGraph :: random_walk(long long start, long long steps) {
    if (start < 0 || start >= num_nodes)
        throw std::out_of_range("Invalid start vertex");

    std::vector<long long> walk;
    walk.reserve(steps + 1);
    walk.push_back(start);

    long long current = start;
    for (long long i = 0; i < steps; i++) {
        if (degree(current) == 0)
        {
            std::cout<<"The seed vertex is an isolated vertex..."<<endl; // isolated vertex, can't continue
            break; 
        }
        current = random_neighbor(current);
        walk.push_back(current);
    }

    return walk;
}

std::vector<double> CSRGraph :: clustering_coeff_sparsification(double p, vector<long long> &walk) {

    long long walk_length = (long long)floor(numEdges() * p);

    long long max_degree = 0;

    for (long long i = 0; i < num_nodes; i++) {
        max_degree = std::max(max_degree, degree(i));
    }

    std::vector<double> triangle_sum(max_degree + 1, 0.0);  // triangle counters per degree
    std::vector<double> node_count(max_degree + 1, 0.0);
    std::vector<double> sparsified_clustering_coeff(max_degree + 1, 0.0);

    // ------------------------------------------------
    // Step 1: Launch a random walk of length walk_length
    // Step 2: Collect distinct edges E and vertices in G_E
    // ------------------------------------------------
    std::set<std::pair<long long, long long>> edge_set; // distinct edges in G_E
    std::unordered_set<long long> vertex_set;           // distinct vertices in G_E

    //std::vector<long long> walk = random_walk(start, walk_length);
    double node_const = (double)(2 * num_edges)/(double)(walk.size());

    for (long long i = 0; i + 1 < (long long)walk.size(); i++) {
        long long u = walk[i];
        long long v = walk[i + 1];

        node_count[degree(walk[i])] += node_const/((double)degree(walk[i]));

        // Canonical form to avoid duplicates
        edge_set.insert({std::min(u, v), std::max(u, v)});
        vertex_set.insert(u);
        vertex_set.insert(v);
    }

    // ------------------------------------------------
    // Step 3: For every vertex v in G_E, count triangles
    //         incident on v using has_edge on G
    //         A triangle {u,v,w} exists in G_E iff all 3
    //         edges {u,v}, {v,w}, {u,w} are in edge_set
    // ------------------------------------------------
    std::vector<double> M(numVertices(), 0.0);

    for (auto& [u, w] : edge_set) {
        // For each edge (u,w) in G_E, find common neighbors v
        // such that (u,v) and (v,w) are also in G_E
        for (long long v : vertex_set) {
            if (v == u || v == w) continue;

            long long lo_uv = std::min(u, v), hi_uv = std::max(u, v);
            long long lo_vw = std::min(v, w), hi_vw = std::max(v, w);

            if (edge_set.count({lo_uv, hi_uv}) &&           //Used to check if {,} is an edge?
                edge_set.count({lo_vw, hi_vw})) {
                // Triangle {u, v, w} found — credit all three vertices
                // Each triangle found 3 times (once per edge), so add 1/3
                M[u] += 1.0 / 3.0;
                M[v] += 1.0 / 3.0;
                M[w] += 1.0 / 3.0;
            }
        }
    }

    // ------------------------------------------------
    // Step 4: Rescale M[v] by 1/p^3 for each vertex v
    // ------------------------------------------------
    double scale = 1.0 / pow(p, 3);
    for (long long v = 0; v < numVertices(); v++)
        M[v] *= scale;


    for(auto v = 0; v < M.size(); v++)
        triangle_sum[degree(v)] += M[v];

    for(auto d = 2; d < triangle_sum.size(); d++)
        sparsified_clustering_coeff[d] += (2.0 * triangle_sum[d]/((d * (d - 1)) * node_count[d]));

    return sparsified_clustering_coeff;
}


std::vector<double> CSRGraph::edge_estimated_degree_clustering_coefficient(double r_frac_edges, double l_frac_edges)
{
    long long r = (long long)(r_frac_edges * num_edges);
    long long l = (long long)(l_frac_edges * num_edges);
    long long max_degree = 0;
    double dR = 0.0;

    std::cout<<"r: "<<r<<", l: "<<l<<std::endl;

    for (long long i = 0; i < num_nodes; i++) {
        max_degree = std::max(max_degree, degree(i));
    }

    std::vector<double> triangle_sum(max_degree + 1, 0.0);  // triangle counters per degree

    std::vector<long long> vtx_degrees, edge_degrees;
    std::vector<std::pair<long long, long long> > edge_samples;

    for(long long i = 0; i < num_nodes; i++)
        vtx_degrees.push_back(degree(i));

    std::discrete_distribution<> vtx_dist(vtx_degrees.begin(), vtx_degrees.end());

    for(long long i = 0; i < r; i++)
    {
        //Sample a vertex proportional to its degree
        long long first_end = vtx_dist(rng);
        //std::cout<<"Vtx: "<<first_end<<std::endl;

        //Sample a uniform random neighbor of this vertex
        long long second_end = random_neighbor(first_end);

        long long edge_degree = (degree(first_end) <= degree(second_end)) ? degree(first_end) : degree(second_end);

        //Add the edge degrees
        dR += edge_degree;

        edge_degrees.push_back(edge_degree);
        edge_samples.push_back(std::make_pair(first_end, second_end));
    }

    std::discrete_distribution<> edge_dist(edge_degrees.begin(), edge_degrees.end());

    for(long long i = 0; i < l; i++)
    {
        //Sample an edge proportional to its degree
        long long sampled_edge_idx = edge_dist(rng);

        long long x = edge_samples[sampled_edge_idx].first;
        long long y = edge_samples[sampled_edge_idx].second;

        //Determine the end point with smaller degree
        long long u = (degree(x) <= degree(y)) ? x : y;
        long long v = (u == x) ? y : x;

        //Sample a neighbor w of u uniformly at random
        long long w = random_neighbor(u);
        
        //Check if w forms a trianggle with v
        if (w != v && has_edge(w, v)) {
            triangle_sum[degree(x)] += 1.0;
            triangle_sum[degree(y)] += 1.0;
        }
    }

    //Rescale the triangle counts
    double scaling_factor = ((((double)num_edges)*(dR))/((2.0 * ((double)l))*(double)r));
    std::cout<<"Scaling factor: "<<scaling_factor<<std::endl;
    for(long long i = 0; i < triangle_sum.size(); i++){
        //std::cout<<"Triangle Count: "<<triangle_sum[i]<<std::endl;
        triangle_sum[i] = (scaling_factor)*((double)triangle_sum[i]/2.0);
    }
        

    return triangle_sum;

}

std::vector<double> CSRGraph :: rw_vertex_estimated_degree_clustering_coefficient(double r_frac_edges, double l_frac_edges, vector<long long> &walk_vertices)
{
    long long r = (long long)(r_frac_edges * num_edges);
    long long l = (long long)(l_frac_edges * num_edges);

    // Arrays indexed by degree
    long long max_degree = 0;

    std::cout<<"Len of RW: "<<r<<", Num Iterations: "<<l<<std::endl;
    for (long long i = 0; i < num_nodes; i++) {
        max_degree = std::max(max_degree, degree(i));
    }

    std::vector<double> triangle_sum(max_degree + 1, 0.0);  // triangle counters per degree
    std::vector<double> node_count(max_degree + 1, 0.0);    // node count per degree


    //Perform a sufficiently long random walk in G
    //vector<long long> walk_vertices = random_walk(seed, r);

    //Rescaling constants to obtain unbiased estimators
    double node_const = (((double)(2.0 * num_edges))/(double) l);
    double triangle_const = (((double) num_edges)/(double) l);

    std::uniform_int_distribution<long long> walk_vertices_dist(0, walk_vertices.size() - 1);
 
    for (long long iter = 0; iter < l; iter++) {

        //Sample a uniform random vertex x on the walk
        long long x = walk_vertices[walk_vertices_dist(rng)];

        //If the sampled vertex is of degree d, then increment node_count[d] by adding 1/d
        node_count[degree(x)] += ((node_const)/(double)degree(x));

        //Sample a uniform random neighbor y of x → edge {x, y}
        long long y = random_neighbor(x);

       // std::cout<<"x: "<<x<<", degree(x): "<<degree(x)<<", y: "<<y<<", degree(y): "<<degree(y)<<std::endl;

        // Determine low and high degree endpoints
        long long u = (degree(x) <= degree(y)) ? x : y;
        long long v = (u == x) ? y : x;

        // Sample a neighbor w of u uniformly at random
        long long w = random_neighbor(u);
        

        //If w is adjacent to v, increment triangle_sum at degree(x) 
        if (has_edge(w, v)) {
            if((((u != v) && (v != w)) &&  u != w))
            {
                triangle_sum[degree(x)] += (triangle_const*((double)degree(u)));    //Need to add degree(u) to the weight because we sampled w uniformly from u's neighbors, so we need to multiply by degree(u) to get the correct inverse probability weight for sampling w.
            }
        }
    }

    // Compute clustering coefficient per degree
    // clustering[d] = triangle_sum[d] / (node_count[d] * (d choose 2))
    std::vector<double> clustering(max_degree + 1, 0.0);
    for (long long d = 2; d <= max_degree; d++) {
        if (node_count[d] == 0.0) continue;

        double d_choose_2 = (double)d * (d - 1) / 2.0;
        //clustering[d] = triangle_sum[d] / (node_count[d] * d_choose_2);
        clustering[d] = triangle_sum[d]/(node_count[d] * d_choose_2);
    }

    return clustering;
}

/*Set scaling factor*/
std::vector<double> CSRGraph::rw_degree_cycle_estimation(double l_frac_edges, vector<long long> &red_nodes, vector<long long> &walk_vertices)
{
    long long max_degree, curr_edge_sample_idx;
    long long l = (long long)(l_frac_edges * num_edges);
    vector<double> walk_edge_weights;
    double total_weight = 0.0;

    for(auto i = 1; i < walk_vertices.size(); i++)
    {
        walk_edge_weights.push_back((double)(degree(walk_vertices[i-1]) * degree(walk_vertices[i])));
        total_weight += ((double)(degree(walk_vertices[i-1]) * degree(walk_vertices[i])));
    }
        

    for (long long i = 0; i < num_nodes; i++) {
        max_degree = std::max(max_degree, degree(i));
    }

    std::vector<double> cycle_sum(max_degree + 1, 0.0);

    std::discrete_distribution<> edge_dist(walk_edge_weights.begin(), walk_edge_weights.end());

    for(auto i = 0; i < l; i++)
    {
        curr_edge_sample_idx = edge_dist(rng);
        long long u = walk_vertices[curr_edge_sample_idx];
        long long v = walk_vertices[curr_edge_sample_idx + 1];

        long long x = random_neighbor(u);
        long long y = random_neighbor(v);


        if((((((u != x) && (u != v)) && ((v != y) && (x != y)))) && (has_edge(x,y))))
        {
            std::map<long long, long long> degree_counts_on_cycle;
            std::vector<long long> red_nodes_on_cycle;

            auto iter = find(red_nodes.begin(), red_nodes.end(), u);
            if(iter != red_nodes.end())
            {
                red_nodes_on_cycle.push_back(u);
                red_nodes_on_cycle.push_back(y);
            }
            else
            {
                red_nodes_on_cycle.push_back(v);
                red_nodes_on_cycle.push_back(x);
            }

            for(auto i = 0; i < red_nodes_on_cycle.size(); i++)
            {
                auto it = degree_counts_on_cycle.find(degree(red_nodes_on_cycle[i]));

                if(it == degree_counts_on_cycle.end())
                    degree_counts_on_cycle[degree(red_nodes_on_cycle[i])] = 0;

                degree_counts_on_cycle[degree(red_nodes_on_cycle[i])] += 1;
                
            }

            for(auto itr = degree_counts_on_cycle.begin(); itr != degree_counts_on_cycle.end(); ++itr) cycle_sum[itr->first] += (double)(1.0/(itr->second));

        }
    }
    double scaling_factor = (double)((num_edges * total_weight)/(l * walk_vertices.size())); //Set this scaling factor properly if needed

    for(auto i = 0; i < cycle_sum.size(); i++)
        cycle_sum[i] = (long long)(scaling_factor * cycle_sum[i]);

    return cycle_sum;
}

//Figure c1, c2, c3
std::vector<double> CSRGraph::estimated_degree_clustering_coefficient(double frac_nodes, double frac_edges) {

    long long r = (long long)(frac_nodes * num_nodes);
    long long l = (long long)(frac_edges * num_edges);

    // Arrays indexed by degree
    long long max_degree = 0;
    std::cout<<"r: "<<r<<", l: "<<l<<std::endl;
    for (long long i = 0; i < num_nodes; i++) {
        max_degree = std::max(max_degree, degree(i));
    }

    std::vector<double> triangle_sum(max_degree + 1, 0.0);  // triangle counters per degree
    std::vector<double> node_count(max_degree + 1, 0.0);    // node count per degree

    // Step 1: Sample r nodes uniformly at random (with replacement) → multiset R
    std::uniform_int_distribution<long long> node_dist(0, num_nodes - 1);
    std::vector<long long> R(r);
    for (long long i = 0; i < r; i++) {
        R[i] = node_dist(rng);
    }

    // Step 2: Compute total degree dR of sampled nodes
    long long dR = 0;
    for (long long x : R) {
        dR += degree(x);
    }

    if (dR == 0) {
        return std::vector<double>(max_degree + 1, 0.0);
    }

    // Step 3: For each node x in R, add c2 to node_count at position degree(x)
    // c2 = num_nodes / r  (inverse sampling probability for nodes)
    double c2 = (double)num_nodes / (double)r;
    for (long long x : R) {
        node_count[degree(x)] += c2;
       // std::cout<<"Degree: "<<degree(x)<<std::endl;
    }

    // Steps 4-7: Repeat l times
    // c1 = (num_nodes * dR) / (r * l)  (inverse sampling probability for edges)
    double c1 = ((double)num_nodes * (double)dR) / ((double)r * (double)l);

    for (long long iter = 0; iter < l; iter++) {

        // Step 4: Sample a node x from R with probability proportional to its degree
        std::uniform_int_distribution<long long> degree_dist(0, dR - 1);
        long long dart = degree_dist(rng);
        long long cumulative = 0;
        long long x = R[0];
        for (long long node : R) {
            cumulative += degree(node);
            if (dart < cumulative) {
                x = node;
                break;
            }
        }

        // Step 5: Sample a uniform random neighbor y of x → edge {x, y}
        long long y = random_neighbor(x);

       // std::cout<<"x: "<<x<<", degree(x): "<<degree(x)<<", y: "<<y<<", degree(y): "<<degree(y)<<std::endl;

        // Determine low and high degree endpoints
        long long u = (degree(x) <= degree(y)) ? x : y;
        long long v = (u == x) ? y : x;

        // Step 6: Sample a neighbor w of u uniformly at random
        long long w = random_neighbor(u);
        

        // Step 7: If w is adjacent to v, increment triangle_sum at degree(x)
        if (w != v && has_edge(w, v)) {
            triangle_sum[degree(x)] += (c1*((double)degree(u)));    //Need to add degree(u) to the weight because we sampled w uniformly from u's neighbors, so we need to multiply by degree(u) to get the correct inverse probability weight for sampling w.
        }
    }

    // Step 8: Compute clustering coefficient per degree
    // clustering[d] = triangle_sum[d] / (node_count[d] * (d choose 2))
    std::vector<double> clustering(max_degree + 1, 0.0);
    for (long long d = 2; d <= max_degree; d++) {
        if (node_count[d] == 0.0) continue;

        double d_choose_2 = (double)d * (d - 1) / 2.0;
        //clustering[d] = triangle_sum[d] / (node_count[d] * d_choose_2);
        clustering[d] = triangle_sum[d];
    }

    return clustering;
}

std::vector<double> CSRGraph:: degree_4cycle_distribution(vector<long long> &red_vertices) {
    // Find max degree
    long long max_deg = 0;
    for (long long v = 0; v < num_nodes; v++) {
        max_deg = std::max(max_deg, degree(v));
    }

    std::vector<double> C(max_deg + 1, 0.0);

    for (long long v = 0; v < num_nodes; v++) {
        long long deg_v = degree(v);
        const long long* nb_v_begin = neighbors_begin(v);
        const long long* nb_v_end   = neighbors_end(v);

        for (const long long* px = nb_v_begin; px != nb_v_end; ++px) {
            long long x = *px;

            for (const long long* pz = px + 1; pz != nb_v_end; ++pz) {
                long long z = *pz;

                const long long* ix = neighbors_begin(x);
                const long long* ex = neighbors_end(x);
                const long long* iz = neighbors_begin(z);
                const long long* ez = neighbors_end(z);

                while (ix != ex && iz != ez) {
                    if (*ix < *iz) {
                        ++ix;
                    } else if (*iz < *ix) {
                        ++iz;
                    } else {
                        long long y = *ix;

                        if (y != v) {
                            long long c = 0;


                            auto it = find(red_vertices.begin(), red_vertices.end(), v);

                            if(it != red_vertices.end())
                            {
                                if (deg_v == deg_v) c++;          // v, always true
                                //if (degree(x) == deg_v) c++;
                                if (degree(y) == deg_v) c++;
                                //if (degree(z) == deg_v) c++;
                            }

                            C[deg_v] += 1.0 / (double)(c);
                        }

                        ++ix;
                        ++iz;
                    }
                }
            }
        }
    }

    return C;
}

/*
void CSRGraph::compare_clustering_coefficients(const std::vector<double>& exact,
    const std::vector<double>& estimated)
{
    long long max_degree = (long long)exact.size() - 1;

    // Step 1: Build bin endpoints as powers of 1.5
    std::vector<long long> bin_endpoints;
    double boundary = 1.0;
    while ((long long)boundary <= max_degree) {
        bin_endpoints.push_back((long long)boundary);
        boundary *= 1.5;
    }
    bin_endpoints.push_back(max_degree + 1); // sentinel for last bin

    long long num_bins = (long long)bin_endpoints.size() - 1;

    // Step 2: For each bin, sum clustering coefficients and count valid degrees
    std::vector<double> exact_bin_avg(num_bins, 0.0);
    std::vector<double> estimated_bin_avg(num_bins, 0.0);
    std::vector<long long> bin_count(num_bins, 0);

    for (long long b = 0; b < num_bins; b++) {
        long long lo = bin_endpoints[b];
        long long hi = bin_endpoints[b + 1]; // exclusive

        for (long long d = lo; d < hi && d <= max_degree; d++) {
            if (d < 2) continue;
            exact_bin_avg[b]     += exact[d];
            estimated_bin_avg[b] += estimated[d];
            bin_count[b]++;
        }
        cout<<"Bin Count: "<<bin_count[b]<<endl;

        if (bin_count[b] > 0) {
            exact_bin_avg[b]     /= bin_count[b];
            estimated_bin_avg[b] /= bin_count[b];
        }
    }

    // Step 3: For every bin, compute relative error and print
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Bin\t\tDegree Range\t\tExact\t\tEstimated\tRelative Error\n";
    std::cout << std::string(90, '-') << "\n";

    std::vector<double> relative_errors; // collect for median/mean

    for (long long b = 0; b < num_bins; b++) {
        if (bin_count[b] == 0) continue;

        long long lo = bin_endpoints[b];
        long long hi = bin_endpoints[b + 1] - 1;

        double relative_error = 0.0;
        if (exact_bin_avg[b] != 0.0) {
            relative_error = std::abs(exact_bin_avg[b] - estimated_bin_avg[b])
                             / exact_bin_avg[b];
        }

        relative_errors.push_back(relative_error);

        std::cout << "Bin " << b
                  << "\t\t[" << lo << ", " << hi << "]"
                  << "\t\t" << exact_bin_avg[b]
                  << "\t\t" << estimated_bin_avg[b]
                  << "\t\t" << relative_error << "\n";
    }

    // Step 4: Compute and print mean and median of per-bin relative errors
    if (relative_errors.empty()) {
        std::cout << "\nNo valid bins to summarize.\n";
        return;
    }

    // Mean
    double mean = 0.0;
    for (double e : relative_errors) mean += e;
    mean /= (double)relative_errors.size();

    // Median
    std::vector<double> sorted_errors = relative_errors;
    std::sort(sorted_errors.begin(), sorted_errors.end());
    double median = 0.0;
    long long n = (long long)sorted_errors.size();
    if (n % 2 == 1) {
        median = sorted_errors[n / 2];
    } else {
        median = (sorted_errors[n / 2 - 1] + sorted_errors[n / 2]) / 2.0;
    }

    std::cout << "\n" << std::string(90, '-') << "\n";
    std::cout << "Summary over " << n << " valid bins:\n";
    std::cout << "  Mean   relative error : " << mean   << "\n";
    std::cout << "  Median relative error : " << median << "\n";
}*/

void CSRGraph::compare_clustering_coefficients(
    const std::vector<double>& exact,
    const std::vector<double>& estimated,
    const std::string& output_file)
{
    long long max_degree = (long long)exact.size() - 1;

    const double alpha = 1.1;  // log base

    // ============================================================
    // Step 1: compute bin index for each degree
    // bin(d) = floor(log(d)/log(alpha))
    // ============================================================

    vector<double> node_counts(max_degree +1, 0.0);
    for(auto i = 0; i < num_nodes; i++)
        node_counts[degree(i)] += 1;

    //for(auto i = 0; i < num_nodes; i++)
    //    std::cout<<"Degree: "<<i<<", count: "<<node_counts[i]<<endl;

    auto get_bin = [&](long long d) -> long long {
        if (d < 2) return -1;
        return (long long)(std::log((double)d) / std::log(alpha));
    };

    long long max_bin = get_bin(max_degree) + 1;

    std::vector<double> exact_bin_sum(max_bin + 1, 0.0);
    std::vector<double> estimated_bin_sum(max_bin + 1, 0.0);
    std::vector<long long> bin_count(max_bin + 1, 0);

    // ============================================================
    // Step 2: accumulate into bins
    // ============================================================

    for (long long d = 2; d <= max_degree; d++) {
        long long b = get_bin(d);
        if (b < 0) continue;
        //cout<<"Node count at degree : "<<d<<" is : "<<node_counts[d]<<endl;
        exact_bin_sum[b]     += exact[d];
        estimated_bin_sum[b] += estimated[d];
        bin_count[b] += node_counts[d];
        //cout<<"Updated val for bin "<<b<<" is: "<<bin_count[b]<<"resulting from degree "<<d<<endl;
    }

    //for(auto i = 0; i < bin_count.size(); i++)
    //    cout<<"Bin: "<<i<<", Count: "<<bin_count[i]<<endl;

    // ============================================================
    // Step 3: print results
    // ============================================================

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Bin\tExactAvg\tEstimatedAvg\tRelativeError\n";
    std::cout << "-----------------------------------------------\n";

    std::vector<double> relative_errors;

    for (long long b = 0; b <= max_bin; b++) {

        if (bin_count[b] == 0) continue;

        double exact_avg = exact_bin_sum[b] / bin_count[b];
        double est_avg   = estimated_bin_sum[b] / bin_count[b];

        double rel_err = 0.0;
        if (exact_avg != 0.0) {
            rel_err = std::abs(exact_avg - est_avg) / exact_avg;
        }

        relative_errors.push_back(rel_err);

        double lo = std::pow(alpha, b);
        double hi = std::pow(alpha, b + 1);

        std::cout << b
                  << "\t[" << lo << "," << hi << "]"
                  << "\t" << exact_avg
                  << "\t" << est_avg
                  << "\t" << rel_err << "\n";
    }

    // ============================================================
    // Step 4: summary stats
    // ============================================================

    if (relative_errors.empty()) return;

    double mean = 0.0;
    for (double x : relative_errors) mean += x;
    mean /= relative_errors.size();

    std::sort(relative_errors.begin(), relative_errors.end());

    double median = relative_errors[relative_errors.size() / 2];

    std::cout << "\nSummary:\n";
    std::cout << "Mean relative error   : " << mean << "\n";
    std::cout << "Median relative error : " << median << "\n";

    std::ofstream out(output_file);

    if (!out.is_open()) {
        std::cerr << "Error: could not open file " << output_file << std::endl;
        return;
    }

    out << std::fixed << std::setprecision(6);

    for (long long b = 0; b <= max_bin; b++) {
        if (bin_count[b] == 0) continue;

        double exact_avg = exact_bin_sum[b] / bin_count[b];
        double est_avg   = estimated_bin_sum[b] / bin_count[b];

        out << exact_avg << " " << est_avg << "\n";
    }

    out.close();
}