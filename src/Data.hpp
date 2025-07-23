#pragma once

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>

class Data {
    private: 
        Eigen::MatrixXd D; // distances matrix

    public: 
        Data(const Eigen::MatrixXd& distances) : D(distances) {}

        // Function to get the distances matrix
        const Eigen::MatrixXd& getDistances() const {
            return D;
        }

        // Function to set the distances matrix
        void setDistances(const Eigen::MatrixXd& distances) {
            D = distances;
        }

        void print() const;
};