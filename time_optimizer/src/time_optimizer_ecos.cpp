#include <iostream>
#include <string>
#include "time_optimizer_ecos.h"
#include "mosek.h"

namespace ecos_sol {

MinimumTimeOptimizer::MinimumTimeOptimizer(){};

MinimumTimeOptimizer::~MinimumTimeOptimizer(){};

int MinimumTimeOptimizer::MinimumTimeGeneration( 
    const Trajectory & traj, const double & maxVel, const double & maxAcc, 
    const double & maxdAcc, const double & d_s, const double &rho) {
    /* minimum time physical feasible trajectory time allocator. */
    /* objective is to generate motion as fast as possible within the physical limitaion (vel, acc and jerk). */
    _P          = traj.getP();
    _T          = traj.getT();
    _poly_num1D = traj.getO();

    uint n_segments  = _P.rows();
    std::cout << "n_segments: " << n_segments << std::endl;

    // int num_a[n_segments], num_b[n_segments], num_c[n_segments], num_d[n_segments];

    std::vector<Eigen::VectorXd> s_list;
    std::vector<uint> k_list;
    uint num_X = 0;
    uint sum_K = 0;
    for(int k = 0; k < n_segments; k++) {
        double duration  = _T(k);
        uint K = ceil(duration / d_s);

        sum_K = sum_K + K;
        k_list.push_back(K);
        Eigen::VectorXd s_k(K + 1);
        for(int i = 0; i < K + 1; i++) {
            s_k(i) = i * d_s;
        }

        s_list.push_back(s_k);

        // num_a[k] = K;
        // num_b[k] = K + 1;
        // num_c[k] = K + 1;
        // num_d[k] = K;
        // num_X = num_X + num_a[k] + num_b[k] + num_c[k] + num_d[k];
        // std::cout << "na: " << num_a[k] << 
        //            "\tnb: " << num_b[k] << 
        //            "\tnc: " << num_c[k] << 
        //            "\tnd: " << num_d[k] << 
        //            "\tnk: " << K << std::endl;
    }

    // Declare the class that helps in finding the X indexes
    uint m = n_segments;
    ecos_sol::variable_set var_set(0, m, k_list);
    uint n_variables = var_set.final_index_ + 1;
    uint size_a = var_set.a_.final_index_ - var_set.a_.initial_index_ + 1;
    uint size_b = var_set.b_.final_index_ - var_set.b_.initial_index_ + 1;
    uint size_c = var_set.c_.final_index_ - var_set.c_.initial_index_ + 1;
    uint size_d = var_set.d_.final_index_ - var_set.d_.initial_index_ + 1;
    // std::cout << size_a << " " << size_b << " " << size_c << " " << size_d << std::endl;

    // Get the number of inequality constraints
    uint n_rot_cones = sum_K;
    uint n_b2c = sum_K + (m + 1);
    uint n_b_positive = sum_K + (m + 1);
    uint n_vel_bounds = 6*sum_K + 6*(m + 1);
    uint n_acc_bounds = 6*sum_K;
    uint n_acc_rate_limiter = 6*sum_K - 6;
    uint n_bcs_acc = 12;

    // Number of equality constraints
    uint n_vel_continuity = 3*m - 1;
    uint n_bcs_vel = 6;

    // Get total number of equalities/inequalities
    uint n_ineq = n_rot_cones + n_b2c + n_b_positive + n_vel_bounds +
                  n_acc_bounds + n_acc_rate_limiter + n_bcs_acc;
    uint n_eq   = n_vel_continuity + n_bcs_vel;
    std::cout << "Number of inequalities: " << n_ineq << std::endl;
    std::cout << "Number of equalities: "   <<   n_eq << std::endl;
    std::cout << "Number of variables: "    << n_variables << std::endl;

    Eigen::MatrixXd G(n_ineq, n_variables);
    Eigen::VectorXd h(n_ineq);
    Eigen::MatrixXd A(n_eq, n_variables);
    Eigen::VectorXd b(n_eq);

    // Current index in the matrix G
    uint cur_index = 0;
    uint n_orthants = 0;


    // Start by setting positive orthants ------------------------------

    // Setup the rotated cone constraints

    // Setup the relation between b and c

    // Setup positiveness of b
    for (uint i = 0; i < n_segments; i++) {
        uint max_K = var_set.b_.segments[i].final_index_ - var_set.b_.segments[i].initial_index_;
        for (uint j = 0; j <= max_K; j++) {
            uint bk_index = var_set.get_index(var_names::b, i, j);
            G(cur_index, bk_index) = 1;
            h(cur_index) = 0;
            cur_index++;
            n_orthants++;
            // std::cout << var_set.get_index(var_names::b, i, j) << std::endl;
        }
    }



    return 1;
}

}  // namespace scs_sol