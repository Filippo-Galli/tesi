#pragma once

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <string>

class Data {
    private:     
        Eigen::MatrixXd D; // distances matrix
        int n; // number of points

        Eigen::VectorXi allocations; // cluster allocations of each point
        int K; // number of clusters
        Eigen::VectorXi cluster_sizes; // sizes of each cluster

    public: 
        Data(const Eigen::MatrixXd& distances, const std::string & first_allocation = "all-in-one");

        // Getters
        double get_distance(int i, int j) const;
        int get_n() const { return n; }
        int get_K() const { return K; }
        const Eigen::VectorXi& get_allocations() const { return allocations; }
        double get_cluster_size(unsigned cluster_index) const {
            return cluster_index < K ? cluster_sizes(cluster_index) : 0;
        }
        double get_cluster_assignment(int index) const { return allocations(index); }
        
        Eigen::VectorXi get_cluster_assignments(int cluster) const;
        
        // Setters
        void set_allocation(int index, int cluster);

        void update_cluster_sizes();
        void update_cluster_sizes(unsigned cluster);
        
        void print() const;
};