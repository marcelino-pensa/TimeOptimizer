#ifndef P4_MAKE_SPARSE_H_
#define P4_MAKE_SPARSE_H_

#include <algorithm> 
#include <iostream> 
#include <vector>
#include <Eigen/Dense>

// Code obtained from:
// https://www.geeksforgeeks.org/sparse-matrix-representations-set-3-csr/

typedef std::vector<std::vector<int> > matrix;
typedef std::vector<int> vi;

// Utility Function to print a Matrix 
void printMatrix(const std::vector<std::vector<int> >& M);

// Utility Function to print A, IA, JA vectors 
// with some decoration. 
void printVector(const std::vector<int> &V, char* msg);

void printVector(const std::vector<double> &V, char* msg);
  
// Generate the three vectors A, IA, JA  
void make_sparse_crs(const std::vector<std::vector<int> >& M);

// Generate the three vectors A, IA, JA  
void make_sparse_ccs(const std::vector<std::vector<int> >& M);

void make_sparse_crs(const Eigen::MatrixXd &M, std::vector<double> *A,
                     std::vector<int> *IA, std::vector<int> *JA);

// Generate the three vectors A, IA, JA  
void make_sparse_ccs(const Eigen::MatrixXd &M, std::vector<double> *A,
                     std::vector<int> *IA, std::vector<int> *JA);

#endif  // P4_MAKE_SPARSE_H_