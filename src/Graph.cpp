#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <random>
#include <stdexcept>
#include "Graph.h"


//Need to update this function
std::vector<double> CSRGraph::actual_degree_clustering_coefficient() {
    
    // Arrays of size num_nodes
    std::vector<long long> triangle_count(num_nodes, 0);  // counter for each vertex
    std::vector<long long> node_degree(num_nodes, 0);     // degree for each vertex

    // Fill degree array
    for (long long i = 0; i < num_nodes; i++) {
        node_degree[i] = degree(i);
    }

    // For each edge (u, v)
    for (long long u = 0; u < num_nodes; u++) {
        for (long long v : neighbors(u)) {
            if (v <= u) continue; // process each edge once

            // Step 1: Determine low degree end
            long long low  = (node_degree[u] <= node_degree[v]) ? u : v;
            long long high = (low == u) ? v : u;
            long long lower_idx = (v < u) ? v : u;

            // Step 2: For every neighbor w of low-degree end,
            // check if high and w are adjacent
            for (long long w : neighbors(low)) {
                if (w == high) continue; // skip the other endpoint
                if (has_edge(high, w)) {
                    if(w > lower_idx) //Added this condition, cross vertify
                    {
                        triangle_count[low]++;
                        triangle_count[high]++;
                        triangle_count[w]++;
                    }
                }
            }
        }
    }

    // Find max degree for sizing
    long long max_degree = *std::max_element(node_degree.begin(), node_degree.end());

    // Step 3: For each degree, count number of nodes of that degree
    std::vector<long long> degree_node_count(max_degree + 1, 0);
    for (long long i = 0; i < num_nodes; i++) {
        degree_node_count[node_degree[i]]++;
    }

    /*for(auto i = 0; i < degree_node_count.size(); i++){
        std::cout<<"Degree: "<<i<<", count: "<<degree_node_count[i]<<std::endl;
    }*/

    // Step 4: For each degree, sum triangle counts of all nodes of that degree
    std::vector<double> degree_triangle_sum(max_degree + 1, 0);
    for (long long i = 0; i < num_nodes; i++) {
        degree_triangle_sum[node_degree[i]] += ((double)triangle_count[i])/2.0;
    }

    // Step 5: For every degree, compute
    // (triangle count for that degree) / (number of nodes of that degree * (degree choose 2))
    std::vector<double> clustering(max_degree + 1, 0.0);
    for (long long d = 2; d <= max_degree; d++) {  // d < 2 means (d choose 2) = 0, skip
        if (degree_node_count[d] == 0) continue;

        double d_choose_2 = (double)d * (d - 1) / 2.0;
        /*clustering[d] = (double)degree_triangle_sum[d] /
                        ((double)degree_node_count[d] * d_choose_2);*/
       clustering[d] = (double)degree_triangle_sum[d];
    }

    return clustering;
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
        triangle_sum[i] = (scaling_factor)*((double)triangle_sum[i]);
    }
        

    return triangle_sum;

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

void CSRGraph::compare_clustering_coefficients(
    const std::vector<double>& exact,
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
}