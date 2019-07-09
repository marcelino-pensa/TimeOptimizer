#ifndef _TIME_OPTIMIZER_SCS_H_
#define _TIME_OPTIMIZER_SCS_H_

#include <Eigen/Dense>
#include <vector>
#include "timeAllocator.h"
#include "trajectory_base.h"
#include "make_sparse.h"
#include "ros/ros.h"
#include "ecos.h"

// Erase libraries
#include <iostream>
#include <fstream>

namespace ecos_sol {

enum class var_names {a, b, c, d};

enum class var_direction {x, y, z};

// Segment goes from initial_index_ to initial_index_+K_
class segment {
 public:
    uint initial_index_;
    uint K_;
    uint final_index_;
    segment() {}
    segment(const uint &initial_index, const uint &K) {
        initial_index_ = initial_index;
        K_ = K;
        final_index_ = initial_index + K;
    }
};

// The solution variables are all stacked in a vector X = [a' b' c' d']'
// The classes below helps returning the index of an entry within X based on
// segment number (0 to m), index number (0 to K), and variable (a, b, c, or d)
class variable {
public: 
    uint initial_index_;
    uint m_;
    uint n_segments_;
    uint final_index_;
    std::vector<segment> segments;

    variable() {}
    variable (const uint &initial_index, const uint &m,
              const std::vector<uint> &K) {
        initial_index_ = initial_index;
        m_ = m;
        n_segments_ = m + 1;
        uint current_index = initial_index;
        for (uint i = 0; i < K.size(); i++) {
            segments.push_back(segment(current_index, K[i]));
            current_index = current_index + K[i] + 1;
        }
        final_index_ = current_index - 1;
    }

    // segment_index - from 0 to m_
    // entry_index - from 0 to K_
    uint get_index(const uint &segment_index, const uint &entry_index) {
        if (segment_index > m_) {
            std::cout << "get_index error: segment_index is above maximum bound" << std::endl;
            return 0;
        }
        if (entry_index > segments[segment_index].K_) {
            std::cout << "get_index error: entry_index is above maximum bound" << std::endl;
            return 0;
        }

        const uint ret = segments[segment_index].initial_index_ + entry_index;
        if (ret < segments[segment_index].initial_index_) {
            std::cout << "get_index error: index is below initial_index" << std::endl;
            return 0;
        } else if (ret > segments[segment_index].final_index_) {
            std::cout << "get_index error: index is above final_index" << std::endl;
            return 0;
        }

        return ret;
    }

    void index_to_ik(const uint &index, std::string *i_k) const {
        const uint index_from_start = index - initial_index_;
        uint i, k;

        // Find the segment to which the  index belongs to
        uint cur_segment = 0;
        uint last_segment_index = 0;
        last_segment_index = last_segment_index + segments[cur_segment].K_;
        while (index_from_start > last_segment_index) {
            cur_segment++;
            last_segment_index = last_segment_index + segments[cur_segment].K_ + 1;
        }
        i = cur_segment;
        k = index - segments[cur_segment].initial_index_;
        *i_k = std::to_string(i) + "_" + std::to_string(k);
    }
};

class variable_set {
public:
    variable a_, b_, c_, d_;
    uint initial_index_;
    uint final_index_;
    variable_set() {}
    variable_set (const uint &initial_index, const uint &m,
                  const std::vector<uint> &K) {
        std::vector<uint> K_minus_1;
        for (uint i = 0; i < K.size(); i++) {
            K_minus_1.push_back(K[i] - 1);
        }
        initial_index_ = initial_index;
        a_ = variable(initial_index, m, K_minus_1);
        b_ = variable(a_.final_index_+1, m, K);
        c_ = variable(b_.final_index_+1, m, K);
        d_ = variable(c_.final_index_+1, m, K_minus_1);
        final_index_ = d_.final_index_;
    }

    uint get_index (const var_names &var_name, const uint i, const uint k) {
        switch(var_name) {
            case var_names::a:
                return a_.get_index(i, k);
            case var_names::b:
                return b_.get_index(i, k);
            case var_names::c:
                return c_.get_index(i, k);
            case var_names::d:
                return d_.get_index(i, k);
        }
    }

    void index_to_var_ik(const uint &index, std::string *var_ik) const {
        std::string i_k;
        if ((index >= a_.initial_index_)  && ((index <= a_.final_index_))) {
            a_.index_to_ik(index, &i_k);
            *var_ik = "a_" + i_k;
        } else if ((index >= b_.initial_index_)  && ((index <= b_.final_index_))) {
            b_.index_to_ik(index, &i_k);
            *var_ik = "b_" + i_k;
        } else if ((index >= c_.initial_index_)  && ((index <= c_.final_index_))) {
            c_.index_to_ik(index, &i_k);
            *var_ik = "c_" + i_k;
        } else if ((index >= d_.initial_index_)  && ((index <= d_.final_index_))) {
            d_.index_to_ik(index, &i_k);
            *var_ik = "d_" + i_k;
        } else {
            std::cout << "Index out of bounds!" << std::endl;
        }
    }
};

class MinimumTimeOptimizer {
private:
        Eigen::MatrixXd _P;  // recording the polynomial's coefficients for further evaluation
        Eigen::VectorXd _T;  // recording the polynomial's time durations
        // int _seg_num, _poly_num1D;

        // double _objective;
        Allocator * time_allocator; // for return the final result to high-level planer

public:
        MinimumTimeOptimizer();
        ~MinimumTimeOptimizer();

        int MinimumTimeGeneration( const Trajectory & traj,
                                    const double & maxVel, 
                                    const double & maxAcc,
                                    const double & maxdAcc,
                                    const double & d_s,
                                    const double &rho); 
        bool check_exit_flag (const idxint &exitflag);

    void vector_to_string(const std::vector<idxint> &vec, const std::string &type,
                      const std::string &var_name, std::string *out);
    void vector_to_string(const std::vector<double> &vec, const std::string &type,
                      const std::string &var_name, std::string *out);
    void vector_to_string(const Eigen::VectorXd &vec, const std::string &type,
                      const std::string &var_name, std::string *out);
    void save_matrices_to_file (
        const uint &n_variables, const uint &size_G, const uint &n_eq,
        const uint &n_orthants, const uint &n_soc, const std::vector<idxint> &soc_size,
        const std::vector<double> &sG, const std::vector<idxint> &JG, const std::vector<idxint> &IG,
        const std::vector<double> &sA, const std::vector<idxint> &JA, const std::vector<idxint> &IA,
        const Eigen::VectorXd &c, const Eigen::VectorXd &h, const Eigen::VectorXd &b);
    void save_socp_constraints_to_files (
        const Eigen::VectorXd &c, const Eigen::MatrixXd &A,
        const Eigen::VectorXd &b, const Eigen::MatrixXd &G,
        const Eigen::VectorXd &h, const uint &n_orthants,
        const std::vector<idxint> &soc_size,
        const ecos_sol::variable_set &var_set);
    void check_constraints_compliance (
        const uint &n_variables, const uint &n_orthants, 
        const uint &n_eq, const pfloat *sol_X,
        const std::vector<std::vector<pfloat>> &sol_a,
        const std::vector<std::vector<pfloat>> &sol_b,
        const std::vector<std::vector<pfloat>> &sol_c,
        const std::vector<std::vector<pfloat>> &sol_d,
        const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
        const Eigen::MatrixXd &G, const Eigen::VectorXd &h);

        Allocator * GetTimeAllocation() {return time_allocator;}
};

}


#endif  // _TIME_OPTIMIZER_SCS_H_
