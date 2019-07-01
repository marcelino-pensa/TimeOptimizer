#ifndef _TIME_OPTIMIZER_SCS_H_
#define _TIME_OPTIMIZER_SCS_H_

#include <Eigen/Dense>
#include <vector>
#include "timeAllocator.h"
#include "trajectory_base.h"

namespace scs_sol {

enum class var_names {a, b, c, d};

enum class var_direction {x, y, z};

// The solution variables are all stacked in a vector X = [a' b' c' d']'
// The classes below helps returning the index of an entry within X based on
// segment number (0 to m), index number (0 to K), and variable (a, b, c, or d)
class variable {
public: 
    uint initial_index_;
    uint m_;
    uint K_;  // segments go from 0 to K_
    uint n_segments_;
    uint n_entries_per_segment_;
    uint final_index_;

    variable() {}
    variable (const uint &initial_index, const uint &m, const uint &K) {
        initial_index_ = initial_index;
        m_ = m;
        n_segments_ = m + 1;
        K_ = K;
        n_entries_per_segment_ = K + 1;
        final_index_ = initial_index + n_entries_per_segment_*n_segments_ - 1;
    }

    // segment_index - from 0 to m_
    // entry_index - from 0 to K_
    uint get_index(const uint &segment_index, const uint &entry_index) {
        if (segment_index > m_) {
            std::cout << "get_index error: segment_index is above maximum bound" << std::endl;
            return 0;
        }
        if (entry_index > K_) {
            std::cout << "get_index error: entry_index is above maximum bound" << std::endl;
            return 0;
        }

        const uint ret = initial_index_ + entry_index + n_entries_per_segment_*segment_index;
        if (ret < initial_index_) {
            std::cout << "get_index error: index is below initial_index" << std::endl;
            return 0;
        } else if (ret > final_index_) {
            std::cout << "get_index error: index is below initial_index" << std::endl;
            return 0;
        }

        return ret;
    }
};

class variable_set {
public:
    variable a_, b_, c_, d_;
    uint final_index_;
    variable_set() {}
    variable_set (const uint &initial_index, const uint &m, const uint K) {
        a_ = variable(initial_index, m, K-1);
        b_ = variable(a_.final_index_+1, m, K);
        c_ = variable(b_.final_index_+1, m, K);
        d_ = variable(c_.final_index_+1, m, K-1);
        final_index_ = d_.final_index_;
    }

    uint get_index (const var_names &var_name, const uint m, const uint K) {
        if (var_name == var_names::a) {
            return a_.get_index(m, K);
        } else if (var_name == var_names::b) {
            return b_.get_index(m, K);
        } else if (var_name == var_names::c) {
            return c_.get_index(m, K);
        } else if (var_name == var_names::d) {
            return d_.get_index(m, K);
        }
    }
};

class variable_set_xyz {
 public:
    variable_set x_, y_, z_;
    uint final_index_;
    variable_set_xyz (const uint &m, const uint K) {
        x_ = variable_set(0, m, K);
        y_ = variable_set(x_.final_index_+1, m, K);
        z_ = variable_set(y_.final_index_+1, m, K);
        final_index_ = z_.final_index_;
    }

    uint get_index (const var_direction &var_xyz, const var_names &var_name,
                    const uint m, const uint K) {
        if (var_xyz == var_direction::x) {
            return x_.get_index(var_name, m, K);
        } else if (var_xyz == var_direction::y) {
            return y_.get_index(var_name, m, K);
        } else if (var_xyz == var_direction::z) {
            return z_.get_index(var_name, m, K);
        }
    }
};

class MinimumTimeOptimizer {
private:
        Eigen::MatrixXd _P;  // recording the polynomial's coefficients for further evaluation
        Eigen::VectorXd _T;  // recording the polynomial's time durations
        int _seg_num, _poly_num1D;

        double _objective;
        Allocator * time_allocator; // for return the final result to high-level planer

public:
        MinimumTimeOptimizer();
        ~MinimumTimeOptimizer();

        int MinimumTimeGeneration( const Trajectory & traj,
                                    const double & maxVel, 
                                    const double & maxAcc,
                                    const double & maxdAcc,
                                    const double & d_s,
                                    const double & rho); 

        Allocator * GetTimeAllocation() {return time_allocator;}
};

}


#endif  // _TIME_OPTIMIZER_SCS_H_
