#ifndef _TIME_OPTIMIZER_SCS_H_
#define _TIME_OPTIMIZER_SCS_H_

#include <Eigen/Dense>
#include <vector>
#include "timeAllocator.h"
#include "trajectory_base.h"
#include "make_sparse.h"
#include "ros/ros.h"

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
};

// class variable_set_xyz {
//  public:
//     variable_set x_, y_, z_;
//     uint final_index_;
//     uint n_variables_;
//     variable_set_xyz (const uint &m, const std::vector<uint> &K) {
//         x_ = variable_set(0, m, K);
//         y_ = variable_set(x_.final_index_+1, m, K);
//         z_ = variable_set(y_.final_index_+1, m, K);
//         final_index_ = z_.final_index_;
//         n_variables_ = final_index_ + 1;
//     }

//     // Get index for i-th segment, k-th entry on x, y, or z (similar to definitions in paper)
//     uint get_index (const var_direction &var_xyz, const var_names &var_name,
//                     const uint i, const uint k) {
//         switch (var_xyz) {
//             case var_direction::x:
//                 return x_.get_index(var_name, i, k);
//             case var_direction::y:
//                 return y_.get_index(var_name, i, k);
//             case var_direction::z:
//                 return z_.get_index(var_name, i, k);
//         }
//     }
// };

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
                                    const double &rho); 

        Allocator * GetTimeAllocation() {return time_allocator;}
};

}


#endif  // _TIME_OPTIMIZER_SCS_H_
