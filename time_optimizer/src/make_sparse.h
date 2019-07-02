#include <algorithm> 
#include <iostream> 
#include <vector>
#include <Eigen/Dense>

// Code obtained from:
// https://www.geeksforgeeks.org/sparse-matrix-representations-set-3-csr/

typedef vector<vector<int> > matrix;
typedef std::vector<int> vi;

// Utility Function to print a Matrix 
void printMatrix(const vector<vector<int> >& M) {
    int m = M.size(); 
    int n = M[0].size(); 
    for (int i = 0; i < m; i++) { 
        for (int j = 0; j < n; j++)  
            std::cout << M[i][j] << " ";         
        std::cout << endl; 
    } 
} 

// Utility Function to print A, IA, JA vectors 
// with some decoration. 
void printVector(const std::vector<int> &V, char* msg) {
    std::cout << msg << "[ "; 
    for_each(V.begin(), V.end(), [](int a) { 
        std::cout << a << " "; 
    }); 
    std::cout << "]" << std::endl; 
}

void printVector(const std::vector<double> &V, char* msg) {
    std::cout << msg << "[ "; 
    for_each(V.begin(), V.end(), [](double a) { 
        std::cout << a << " "; 
    }); 
    std::cout << "]" << std::endl; 
}
  
// Generate the three vectors A, IA, JA  
void make_sparse_crs(const vector<vector<int> >& M) {
    int m = M.size(); 
    int n = M[0].size(), i, j; 
    std::vector<double> A; 
    std::vector<int> IA = { 0 }; // IA matrix has N+1 rows 
    std::vector<int> JA; 
    int NNZ = 0; 
  
    for (i = 0; i < m; i++) { 
        for (j = 0; j < n; j++) { 
            if (M[i][j] != 0) { 
                A.push_back(M[i][j]); 
                JA.push_back(j); 
  
                // Count Number of Non Zero  
                // Elements in row i 
                NNZ++; 
            } 
        } 
        IA.push_back(NNZ); 
    } 
  
    // printMatrix(M);
    // printVector(A, (char*)"A = ");
    // printVector(IA, (char*)"IA = ");
    // printVector(JA, (char*)"JA = ");
}

// Generate the three vectors A, IA, JA  
void make_sparse_ccs(const vector<vector<int> >& M) {
    int m = M[0].size();  // n columns
    int n = M.size();  // n rows
    int i, j;
    std::vector<double> A; 
    std::vector<int> JA = { 0 }; // IA matrix has N+1 rows 
    std::vector<int> IA; 
    int NNZ = 0; 
  
    for (i = 0; i < m; i++) { 
        for (j = 0; j < n; j++) { 
            if (M[j][i] != 0) { 
                A.push_back(M[j][i]); 
                IA.push_back(j); 
  
                // Count Number of Non Zero  
                // Elements in column j
                NNZ++; 
            } 
        } 
        JA.push_back(NNZ); 
    } 
  
    // printMatrix(M); 
    // printVector(A, (char*)"A = "); 
    // printVector(IA, (char*)"IA = "); 
    // printVector(JA, (char*)"JA = "); 
}

void make_sparse_crs(const Eigen::MatrixXd &M, std::vector<double> *A,
                     std::vector<int> *IA, std::vector<int> *JA) {
    int m = M.rows();
    int n = M.cols();
    IA->push_back(0); // IA matrix has N+1 rows 
    int NNZ = 0; 
  
    for (uint i = 0; i < m; i++) { 
        for (uint j = 0; j < n; j++) { 
            if (M(i,j) != 0) { 
                A->push_back(M(i,j)); 
                JA->push_back(j); 
  
                // Count Number of Non Zero  
                // Elements in row i 
                NNZ++; 
            } 
        } 
        IA->push_back(NNZ); 
    } 
}

// Generate the three vectors A, IA, JA  
void make_sparse_ccs(const Eigen::MatrixXd &M, std::vector<double> *A,
                     std::vector<int> *IA, std::vector<int> *JA) {
    int m = M.cols();  // n columns
    int n = M.rows();  // n rows
    JA->push_back(0); // JA matrix has N+1 rows 
    int NNZ = 0; 
  
    for (uint i = 0; i < m; i++) {
        for (uint j = 0; j < n; j++) {
            if (M(j,i) != 0) { 
                A->push_back(M(j,i)); 
                IA->push_back(j); 
  
                // Count Number of Non Zero  
                // Elements in column j
                NNZ++; 
            } 
        } 
        JA->push_back(NNZ); 
    }
    // printVector(*A, (char*)"A = ");
    // printVector(*IA, (char*)"IA = ");
    // printVector(*JA, (char*)"JA = ");
}