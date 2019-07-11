#ifndef P4_MAKE_SPARSE_H_
#define P4_MAKE_SPARSE_H_

#include <algorithm> 
#include <iostream> 
#include <vector>
#include <Eigen/Dense>
#include "glblopts.h"

namespace sparse {
// Code obtained from:
// https://www.geeksforgeeks.org/sparse-matrix-representations-set-3-csr/

class row_val {
 public:
    idxint row_;
    double value_;
    row_val (const idxint &row, const double &value){
        row_ = row;
        value_ = value;
    }
};

class col_val {
 public:
    idxint col_;
    double value_;
    col_val (const idxint &col, const double &value){
        col_ = col;
        value_ = value;
    }
};

// Entry M[i] has the i-th column
// This sparse representation is made to be easily transformed into
//   CCS form while being added row by row
class sp_matrix {
 public:
    std::vector<std::vector<row_val>> Mat_;
    idxint n_rows_, n_cols_;

    sp_matrix (const idxint &n_cols) {
        n_rows_ = 0;
        n_cols_ = n_cols;
        Mat_.resize(n_cols_);
    }

    void add_row(const std::vector<idxint> &columns, 
                 const std::vector<idxint> &values) {
        for (uint i = 0; i < columns.size(); i++) {
            Mat_[columns[i]].push_back(row_val(n_rows_, values[i]));
        }
        n_rows_++;
    }

    void add_row(const std::vector<col_val> &col_values) {
        for (uint i = 0; i < col_values.size(); i++) {
            Mat_[col_values[i].col_].push_back(row_val(n_rows_, col_values[i].value_));
        }
        n_rows_++;
    }

    void to_ccs(std::vector<double> *A, std::vector<idxint> *IA, 
                std::vector<idxint> *JA) {
        JA->push_back(0); // JA matrix has n_cols+1 entries, with the first one being zero
        int NNZ = 0;      // number of non-zero entries
        for (uint i = 0; i < n_cols_; i++) {
            for (uint j = 0; j < Mat_[i].size(); j++) {
                A->push_back(Mat_[i][j].value_); 
                IA->push_back(Mat_[i][j].row_);
  
                // Count Number of Non Zero elements in column j
                NNZ++;
            } 
            JA->push_back(NNZ); 
        }
    }
};

// typedef std::vector<std::vector<int> > matrix;
// typedef std::vector<int> vi;

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
                     std::vector<idxint> *IA, std::vector<idxint> *JA);

// Generate the three vectors A, IA, JA  
void make_sparse_ccs(const Eigen::MatrixXd &M, std::vector<double> *A,
                     std::vector<idxint> *IA, std::vector<idxint> *JA);

} // sparse

#endif  // P4_MAKE_SPARSE_H_