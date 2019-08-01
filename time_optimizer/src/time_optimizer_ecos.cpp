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
    // _poly_num1D = traj.getO();
    double maxJer_s = maxdAcc * d_s;

    uint n_segments  = _P.rows();
    std::cout << "n_segments: " << n_segments << std::endl;
    std::cout << "segment times: " << _T.transpose() << std::endl;

    std::vector<Eigen::VectorXd> t_list;
    std::vector<uint> k_list;
    uint num_X = 0;
    uint sum_K = 0;
    double new_ds; // d_s might change slightly if duration/d_s is not an integer
    for(int j = 0; j < n_segments; j++) {
        double duration  = _T(j);
        uint K = ceil(duration / d_s);

        // Get lists of times
        sum_K = sum_K + K;
        k_list.push_back(K);
        Eigen::VectorXd t_k(K + 1);
        new_ds = duration/K;
        for(int i = 0; i < K + 1; i++) {
            t_k(i) = i * new_ds;
        }
        t_list.push_back(t_k);
    }

    // Declare the class that helps in finding the X indexes
    uint m = n_segments-1;
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
    uint n_vel_bounds = 3*sum_K + 3*(m + 1);
    uint n_acc_bounds = 6*sum_K;
    uint n_acc_rate_limiter = 6*sum_K - 6;
    uint n_bcs_acc = 12;

    // Number of equality constraints
    uint n_vel_continuity = m;
    uint n_b2a = sum_K;
    uint n_bcs_vel = 2;

    // Get total number of equalities/inequalities
    uint n_ineq = n_rot_cones + n_b2c + n_b_positive + n_vel_bounds +
                  n_acc_bounds + n_acc_rate_limiter + n_bcs_acc;
    uint size_G = 3*n_rot_cones + 3*n_b2c + n_b_positive + n_vel_bounds +
                  n_acc_bounds + n_acc_rate_limiter + n_bcs_acc;
    // uint size_G = 3*n_rot_cones + 3*n_b2c + n_b_positive + n_vel_bounds +
    //               n_acc_bounds + n_bcs_acc;
    uint n_eq   = n_vel_continuity + n_b2a;
    // std::cout << "Number of inequalities: " << n_ineq << std::endl;
    // std::cout << "Number of equalities: "   <<   n_eq << std::endl;
    // std::cout << "Number of variables: "    << n_variables << std::endl;

    Eigen::VectorXd c = Eigen::VectorXd::Zero(n_variables);
    // Eigen::MatrixXd G = Eigen::MatrixXd::Zero(size_G, n_variables);
    Eigen::VectorXd h = Eigen::VectorXd::Zero(size_G);
    // Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n_eq, n_variables);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(n_eq);
    sparse::sp_matrix G_sparse(n_variables);
    sparse::sp_matrix A_sparse(n_variables);

    // Current index in the matrix G
    uint cur_G_index = 0;
    uint cur_A_index = 0;
    uint n_orthants = 0;
    uint n_soc = 0;  // Number of second order cones

    // ------------------------------------------------------------------
    // Set cost function ------------------------------------------------
    // ------------------------------------------------------------------
    const uint d_init_index = var_set.d_.initial_index_;
    const uint d_final_index = var_set.d_.final_index_;
    for (uint i = d_init_index; i <= d_final_index; i++) {
        c[i] = 1;
    }

    // ------------------------------------------------------------------
    // Set positive orthants --------------------------------------------
    // ------------------------------------------------------------------
    // Positive orthants are of the type G*x < h <==> h - G*x > 0

    // Set b > 0  <==> -b < 0
    for (uint i = 0; i <= m; i++) {
        uint max_K = var_set.b_.segments[i].K_;
        // std::cout << "maxk b: " << max_K << std::endl;
        for (uint j = 0; j <= max_K; j++) {
            const uint bk_index = var_set.get_index(var_names::b, i, j);
            std::vector<sparse::col_val> col_values;
            col_values.push_back(sparse::col_val(bk_index, -1));
            G_sparse.add_row(col_values);

            // G(cur_G_index, bk_index) = -1;
            h(cur_G_index) = 0;
            cur_G_index++;
            n_orthants++;
            // std::cout << var_set.get_index(var_names::b, i, j) << std::endl;
        }
    }

    // Set velocity bounds
    // f'^2 . b_i_k < Vmax^2
    const double max_vel_sqr = maxVel*maxVel;
    for (uint i = 0; i <= m; i++) {
        uint max_K = var_set.b_.segments[i].K_;
        Eigen::VectorXd time = t_list[i];
        for (uint k = 0; k <= max_K; k++) {
            const uint bk_index = var_set.get_index(var_names::b, i, k);
            const Eigen::Vector3d f_prime = traj.getVel(i, time(k));

            // Add for x, y, and z directions
            for (uint idx = 0; idx < 3; idx++) {
                std::vector<sparse::col_val> col_values;
                col_values.push_back(sparse::col_val(bk_index, f_prime(idx)*f_prime(idx)));
                G_sparse.add_row(col_values);

                // G(cur_G_index, bk_index) = f_prime(idx)*f_prime(idx);
                h(cur_G_index) = max_vel_sqr;
                cur_G_index++;
                n_orthants++;
            }
        }
    }

    // Set acceleration bounds
    //   f' . a_i_k + f'' . b_i_k < a_max
    // - f' . a_i_k - f'' . b_i_k < a_max
    const double half_ds = 0.5*new_ds;
    for (uint i = 0; i <= m; i++) {
        uint max_Ka = var_set.a_.segments[i].K_;
        // std::cout << "maxk b: " << max_K << std::endl;
        Eigen::VectorXd time = t_list[i];
        for (uint k = 0; k <= max_Ka; k++) {
            const uint ak_index = var_set.get_index(var_names::a, i, k);
            const uint bk_index = var_set.get_index(var_names::b, i, k);
            const uint bkplus1_index = bk_index + 1;
            const double eval_time = time(k) + half_ds;
            const Eigen::Vector3d f_prime = traj.getVel(i, eval_time);
            const Eigen::Vector3d f_prime2 = traj.getAcc(i, eval_time);

            // Add for x, y, and z directions (upper bound)
            for (uint idx = 0; idx < 3; idx++) {
                std::vector<sparse::col_val> col_values;
                col_values.push_back(sparse::col_val(ak_index,      f_prime(idx)));
                col_values.push_back(sparse::col_val(bk_index,      0.5*f_prime2(idx)));
                col_values.push_back(sparse::col_val(bkplus1_index, 0.5*f_prime2(idx)));
                G_sparse.add_row(col_values);

                // G(cur_G_index, ak_index) =          f_prime(idx);
                // G(cur_G_index, bk_index) =      0.5*f_prime2(idx);
                // G(cur_G_index, bkplus1_index) = 0.5*f_prime2(idx);
                h(cur_G_index) = maxAcc;
                cur_G_index++;
                n_orthants++;
            }

            // Add for x, y, and z directions (lower bound)
            for (uint idx = 0; idx < 3; idx++) {
                std::vector<sparse::col_val> col_values;
                col_values.push_back(sparse::col_val(ak_index,     -f_prime(idx)));
                col_values.push_back(sparse::col_val(bk_index,     -0.5*f_prime2(idx)));
                col_values.push_back(sparse::col_val(bkplus1_index,-0.5*f_prime2(idx)));
                G_sparse.add_row(col_values);

                // G(cur_G_index, ak_index) =          -f_prime(idx);
                // G(cur_G_index, bk_index) =      -0.5*f_prime2(idx);
                // G(cur_G_index, bkplus1_index) = -0.5*f_prime2(idx);
                h(cur_G_index) = maxAcc;
                cur_G_index++;
                n_orthants++;
            }
        }
    }

    // Set rate limiter
    const double half_plus1_ds = 1.5*new_ds;
    for (uint i = 0; i <= m; i++) {
        uint max_Ka = var_set.a_.segments[i].K_;
        // std::cout << "maxk b: " << max_K << std::endl;
        Eigen::VectorXd time = t_list[i];
        for (uint k = 0; k < max_Ka; k++) {
            const uint ak_index = var_set.get_index(var_names::a, i, k);
            const uint akplus1_index = ak_index + 1;
            const uint bk_index = var_set.get_index(var_names::b, i, k);
            const uint bkplus1_index = bk_index + 1;
            const uint bkplus2_index = bk_index + 2;
            const double eval_time = time(k) + half_ds;
            const double eval_time2 = time(k) + half_plus1_ds;
            const Eigen::Vector3d f_prime = traj.getVel(i, eval_time);
            const Eigen::Vector3d f_prime2 = traj.getAcc(i, eval_time);
            const Eigen::Vector3d f_plus1_prime = traj.getVel(i, eval_time2);
            const Eigen::Vector3d f_plus1_prime2 = traj.getAcc(i, eval_time2);

            // Add for x, y, and z directions (upper bound)
            for (uint idx = 0; idx < 3; idx++) {
                std::vector<sparse::col_val> col_values;
                col_values.push_back(sparse::col_val(ak_index,       f_prime(idx)));
                col_values.push_back(sparse::col_val(akplus1_index, -f_plus1_prime(idx)));
                col_values.push_back(sparse::col_val(bk_index,       0.5*f_prime2(idx)));
                col_values.push_back(sparse::col_val(bkplus1_index,  0.5*f_prime2(idx) - 0.5*f_plus1_prime2(idx)));
                col_values.push_back(sparse::col_val(bkplus2_index, -0.5*f_plus1_prime2(idx)));
                G_sparse.add_row(col_values);

                // G(cur_G_index, ak_index) =           f_prime(idx);
                // G(cur_G_index, akplus1_index) =     -f_plus1_prime(idx);
                // G(cur_G_index, bk_index) =       0.5*f_prime2(idx);
                // G(cur_G_index, bkplus1_index) =  0.5*f_prime2(idx) - 0.5*f_plus1_prime2(idx);
                // G(cur_G_index, bkplus2_index) = -0.5*f_plus1_prime2(idx);
                h(cur_G_index) = maxJer_s;
                cur_G_index++;
                n_orthants++;
            }

            // Add for x, y, and z directions (lower bound)
            for (uint idx = 0; idx < 3; idx++) {
                std::vector<sparse::col_val> col_values;
                col_values.push_back(sparse::col_val(ak_index,      -f_prime(idx)));
                col_values.push_back(sparse::col_val(akplus1_index,  f_plus1_prime(idx)));
                col_values.push_back(sparse::col_val(bk_index,      -0.5*f_prime2(idx)));
                col_values.push_back(sparse::col_val(bkplus1_index, -0.5*f_prime2(idx) + 0.5*f_plus1_prime2(idx)));
                col_values.push_back(sparse::col_val(bkplus2_index,  0.5*f_plus1_prime2(idx)));
                G_sparse.add_row(col_values);

                // G(cur_G_index, ak_index) =          -f_prime(idx);
                // G(cur_G_index, akplus1_index) =      f_plus1_prime(idx);
                // G(cur_G_index, bk_index) =      -0.5*f_prime2(idx);
                // G(cur_G_index, bkplus1_index) = -0.5*f_prime2(idx) + 0.5*f_plus1_prime2(idx);
                // G(cur_G_index, bkplus2_index) =  0.5*f_plus1_prime2(idx);
                h(cur_G_index) = maxJer_s;
                cur_G_index++;
                n_orthants++;
            }
        }
    }

    // Set rate limited between segments
    for (uint i = 1; i <= m; i++) {
        uint max_Ka = var_set.a_.segments[i-1].K_;
        const Eigen::VectorXd time = t_list[i-1];
        const double eval_time = half_ds;
        const double eval_time2 = time(time.size()-1) - half_ds;
        const Eigen::Vector3d f_prime = traj.getVel(i, eval_time);
        const Eigen::Vector3d f_prime2 = traj.getAcc(i, eval_time);
        const Eigen::Vector3d f_minus1_prime = traj.getVel(i-1, eval_time2);
        const Eigen::Vector3d f_minus1_prime2 = traj.getAcc(i-1, eval_time2);
        
        const uint a_iminus1_kminus1_index = var_set.get_index(var_names::a, i-1, max_Ka);
        const uint a_i_0_index = var_set.get_index(var_names::a, i, 0);
        const uint b_iminus1_kminus1_index = var_set.get_index(var_names::b, i-1, max_Ka);
        const uint b_iminus1_k_index = b_iminus1_kminus1_index + 1;
        const uint b_i_0_index = var_set.get_index(var_names::b, i, 0);
        const uint b_i_1_index = var_set.get_index(var_names::b, i, 1);

        // Add for x, y, and z directions (upper bound)
        for (uint idx = 0; idx < 3; idx++) {
            std::vector<sparse::col_val> col_values;
            col_values.push_back(sparse::col_val(a_iminus1_kminus1_index,-f_minus1_prime(idx)));
            col_values.push_back(sparse::col_val(a_i_0_index,             f_prime(idx)));
            col_values.push_back(sparse::col_val(b_iminus1_kminus1_index,-0.5*f_minus1_prime2(idx)));
            col_values.push_back(sparse::col_val(b_iminus1_k_index,      -0.5*f_minus1_prime2(idx)));
            col_values.push_back(sparse::col_val(b_i_0_index,             0.5*f_prime2(idx)));
            col_values.push_back(sparse::col_val(b_i_1_index,             0.5*f_prime2(idx)));
            G_sparse.add_row(col_values);

            // G(cur_G_index, a_iminus1_kminus1_index) =-f_minus1_prime(idx);
            // G(cur_G_index, a_i_0_index) =             f_prime(idx);
            // G(cur_G_index, b_iminus1_kminus1_index) =-0.5*f_minus1_prime2(idx);
            // G(cur_G_index, b_iminus1_k_index) =      -0.5*f_minus1_prime2(idx);
            // G(cur_G_index, b_i_0_index) =             0.5*f_prime2(idx);
            // G(cur_G_index, b_i_1_index) =             0.5*f_prime2(idx);
            h(cur_G_index) = maxJer_s;
            cur_G_index++;
            n_orthants++;
        }

        // Add for x, y, and z directions (lower bound)
        for (uint idx = 0; idx < 3; idx++) {
            std::vector<sparse::col_val> col_values;
            col_values.push_back(sparse::col_val(a_iminus1_kminus1_index, f_minus1_prime(idx)));
            col_values.push_back(sparse::col_val(a_i_0_index,            -f_prime(idx)));
            col_values.push_back(sparse::col_val(b_iminus1_kminus1_index, 0.5*f_minus1_prime2(idx)));
            col_values.push_back(sparse::col_val(b_iminus1_k_index,       0.5*f_minus1_prime2(idx)));
            col_values.push_back(sparse::col_val(b_i_0_index,            -0.5*f_prime2(idx)));
            col_values.push_back(sparse::col_val(b_i_1_index,            -0.5*f_prime2(idx)));
            G_sparse.add_row(col_values);

            // G(cur_G_index, a_iminus1_kminus1_index) = f_minus1_prime(idx);
            // G(cur_G_index, a_i_0_index) =            -f_prime(idx);
            // G(cur_G_index, b_iminus1_kminus1_index) = 0.5*f_minus1_prime2(idx);
            // G(cur_G_index, b_iminus1_k_index) =       0.5*f_minus1_prime2(idx);
            // G(cur_G_index, b_i_0_index) =            -0.5*f_prime2(idx);
            // G(cur_G_index, b_i_1_index) =            -0.5*f_prime2(idx);
            h(cur_G_index) = maxJer_s;
            cur_G_index++;
            n_orthants++;
        }
    }

    // Acceleration boundary conditions (set to zero) -------------------------
    const uint a0_index = var_set.a_.initial_index_;
    const uint b0_index = var_set.b_.initial_index_;
    const uint b1_index = b0_index + 1;
    const uint last_a_index = var_set.a_.final_index_;
    const uint last_b_index = var_set.b_.final_index_;
    const uint b4last_b_index = last_b_index - 1;
    const Eigen::VectorXd last_time = t_list[t_list.size()-1];
    const double final_time = last_time[last_time.size()-1];
    const Eigen::Vector3d f_prime_init  = traj.getVel(0, new_ds);
    const Eigen::Vector3d f_prime2_init = traj.getAcc(0, new_ds);
    const Eigen::Vector3d f_prime_final  = traj.getVel(m, final_time - new_ds);
    const Eigen::Vector3d f_prime2_final = traj.getAcc(m, final_time - new_ds);
    const double init_acc = 0, final_acc = 0;

    // Initial boundary conditions for acceleration (upper bound)
    for (uint idx = 0; idx < 3; idx++) {
        std::vector<sparse::col_val> col_values;
        col_values.push_back(sparse::col_val(a0_index,     f_prime_init(idx)));
        col_values.push_back(sparse::col_val(b0_index, 0.5*f_prime2_init(idx)));
        col_values.push_back(sparse::col_val(b1_index, 0.5*f_prime2_init(idx)));
        G_sparse.add_row(col_values);

        // G(cur_G_index, a0_index) = f_prime_init(idx);
        // G(cur_G_index, b0_index) = 0.5*f_prime2_init(idx);
        // G(cur_G_index, b1_index) = 0.5*f_prime2_init(idx);
        h(cur_G_index) = 0.01 + init_acc;
        cur_G_index++;
        n_orthants++;
    }

    // Initial boundary conditions for acceleration (lower bound)
    for (uint idx = 0; idx < 3; idx++) {
        std::vector<sparse::col_val> col_values;
        col_values.push_back(sparse::col_val(a0_index,     -f_prime_init(idx)));
        col_values.push_back(sparse::col_val(b0_index, -0.5*f_prime2_init(idx)));
        col_values.push_back(sparse::col_val(b1_index, -0.5*f_prime2_init(idx)));
        G_sparse.add_row(col_values);

        // G(cur_G_index, a0_index) = -f_prime_init(idx);
        // G(cur_G_index, b0_index) = -0.5*f_prime2_init(idx);
        // G(cur_G_index, b1_index) = -0.5*f_prime2_init(idx);
        h(cur_G_index) = 0.01 - init_acc;
        cur_G_index++;
        n_orthants++;
    }

    // Final boundary conditions for acceleration (upper bound)
    for (uint idx = 0; idx < 3; idx++) {
        std::vector<sparse::col_val> col_values;
        col_values.push_back(sparse::col_val(last_a_index,       f_prime_final(idx)));
        col_values.push_back(sparse::col_val(b4last_b_index, 0.5*f_prime2_final(idx)));
        col_values.push_back(sparse::col_val(last_b_index,   0.5*f_prime2_final(idx)));
        G_sparse.add_row(col_values);

        // G(cur_G_index, last_a_index) =   f_prime_final(idx);
        // G(cur_G_index, b4last_b_index) = 0.5*f_prime2_final(idx);
        // G(cur_G_index, last_b_index) =   0.5*f_prime2_final(idx);
        h(cur_G_index) = 0.01 + final_acc;
        cur_G_index++;
        n_orthants++;
    }

    // Final boundary conditions for acceleration (lower bound)
    for (uint idx = 0; idx < 3; idx++) {
        std::vector<sparse::col_val> col_values;
        col_values.push_back(sparse::col_val(last_a_index,       -f_prime_final(idx)));
        col_values.push_back(sparse::col_val(b4last_b_index, -0.5*f_prime2_final(idx)));
        col_values.push_back(sparse::col_val(last_b_index,   -0.5*f_prime2_final(idx)));
        G_sparse.add_row(col_values);

        // G(cur_G_index, last_a_index) =   -f_prime_final(idx);
        // G(cur_G_index, b4last_b_index) = -0.5*f_prime2_final(idx);
        // G(cur_G_index, last_b_index) =   -0.5*f_prime2_final(idx);
        h(cur_G_index) = 0.01 - final_acc;
        cur_G_index++;
        n_orthants++;
    }

    // ------------------------------------------------------------------
    // Set second order cones--------------------------------------------
    // ------------------------------------------------------------------
    std::vector<idxint> soc_size;

    // Setup the rotated cone constraints -----------------------------
    // dk*(c_{k+1} + c_k)
    for (uint i = 0; i <= m; i++) {
        uint max_Kd = var_set.d_.segments[i].K_;
        for (uint k = 0; k <= max_Kd; k++) {
            std::vector<sparse::col_val> col_values0, col_values1, col_values2;
            const uint dk_index = var_set.get_index(var_names::d, i, k);
            const uint ck_index = var_set.get_index(var_names::c, i, k);
            const uint ckplus1_index = ck_index + 1;

            col_values0.push_back(sparse::col_val(ck_index,     -1));
            col_values0.push_back(sparse::col_val(ckplus1_index,-1));
            col_values0.push_back(sparse::col_val(dk_index,     -1));
            // col_values1 has an empty row
            col_values2.push_back(sparse::col_val(ck_index,      1));
            col_values2.push_back(sparse::col_val(ckplus1_index, 1));
            col_values2.push_back(sparse::col_val(dk_index,     -1));
            G_sparse.add_row(col_values0);
            G_sparse.add_row(col_values1);
            G_sparse.add_row(col_values2);

            // G(cur_G_index + 0, ck_index) =     -1;
            // G(cur_G_index + 0, ckplus1_index) =-1;
            // G(cur_G_index + 0, dk_index) =     -1;
            // G(cur_G_index + 2, ck_index) =      1;
            // G(cur_G_index + 2, ckplus1_index) = 1;
            // G(cur_G_index + 2, dk_index) =     -1;
            h(cur_G_index + 0) = 0;
            h(cur_G_index + 1) = 2;
            h(cur_G_index + 2) = 0;
            cur_G_index = cur_G_index + 3;
            soc_size.push_back(3);
            n_soc++;
        }
    }

    // Setup the relation between b and c
    for (uint i = 0; i <= m; i++) {
        uint max_Kb = var_set.b_.segments[i].K_;
        for (uint k = 0; k <= max_Kb; k++) {
            std::vector<sparse::col_val> col_values0, col_values1, col_values2;
            const uint bk_index = var_set.get_index(var_names::b, i, k);
            const uint ck_index = var_set.get_index(var_names::c, i, k);

            col_values0.push_back(sparse::col_val(bk_index, -1));
            col_values0.push_back(sparse::col_val(ck_index,  0));
            col_values1.push_back(sparse::col_val(bk_index, -1));
            col_values1.push_back(sparse::col_val(ck_index,  0));
            col_values2.push_back(sparse::col_val(bk_index,  0));
            col_values2.push_back(sparse::col_val(ck_index, -2));
            G_sparse.add_row(col_values0);
            G_sparse.add_row(col_values1);
            G_sparse.add_row(col_values2);

            // G(cur_G_index + 0, bk_index) =-1;
            // G(cur_G_index + 0, ck_index) = 0;
            // G(cur_G_index + 1, bk_index) =-1;
            // G(cur_G_index + 1, ck_index) = 0;
            // G(cur_G_index + 2, bk_index) = 0;
            // G(cur_G_index + 2, ck_index) =-2;
            h(cur_G_index + 0) =  1;
            h(cur_G_index + 1) = -1;
            h(cur_G_index + 2) =  0;
            cur_G_index = cur_G_index + 3;
            soc_size.push_back(3);
            n_soc++;
        }
    }

    // ------------------------------------------------------------------
    // Set equality constraints -----------------------------------------
    // ------------------------------------------------------------------
    // Velocity continuity between every two segments
    for (uint i = 0; i < m; i++) {
        std::vector<sparse::col_val> col_values;
        const uint max_K = var_set.b_.segments[i].K_;
        const uint b_i_K_index = var_set.get_index(var_names::b, i, max_K);
        const uint b_iplus1_0_index = var_set.get_index(var_names::b, i+1, 0);
        col_values.push_back(sparse::col_val(b_i_K_index,       1));
        col_values.push_back(sparse::col_val(b_iplus1_0_index, -1));
        A_sparse.add_row(col_values);
        // A(cur_A_index, b_i_K_index) = 1;
        // A(cur_A_index, b_iplus1_0_index) = -1;
        b(cur_A_index) = 0;
        cur_A_index++;
    }

    // Relation between b and a
    // bi_k1 - bi_k - 2.ds.ai_k = 0
    for (uint i = 0; i <= m; i++) {
        uint max_Ka = var_set.a_.segments[i].K_;
        for (uint k = 0; k <= max_Ka; k++) {
            std::vector<sparse::col_val> col_values;
            const uint ak_index = var_set.get_index(var_names::a, i, k);
            const uint bk_index = var_set.get_index(var_names::b, i, k);
            const uint bk1_index = bk_index + 1;
            col_values.push_back(sparse::col_val(ak_index, -2*new_ds));
            col_values.push_back(sparse::col_val(bk_index, -1));
            col_values.push_back(sparse::col_val(bk1_index, 1));
            A_sparse.add_row(col_values);
            // A(cur_A_index, ak_index) = -2*new_ds;
            // A(cur_A_index, bk_index) = -1;
            // A(cur_A_index, bk1_index) = 1;
            b(cur_A_index) = 0;
            cur_A_index++;
        }
    }

    // std::cout << "added inequalities: " << n_soc+n_orthants << std::endl;
    // std::cout << "added equalities: "   << cur_A_index      << std::endl;

    ros::Time t0 = ros::Time::now();
    std::vector<double> sG, sA;
    std::vector<idxint> IG, JG, IA, JA;
    G_sparse.to_ccs(&sG, &IG, &JG);
    A_sparse.to_ccs(&sA, &IA, &JA);

    ros::Time t1 = ros::Time::now();
    std::cout << "time to make sparse: " << (t1 - t0).toSec() << std::endl;

    // Velocity boundary conditions (set initial/final velocity to zero)
    // This doesn't seem to be needed
    // {
    //     const uint max_K = var_set.b_.segments[m].K_;
    //     const uint init_b_index  = var_set.get_index(var_names::b, 0, 0);
    //     const uint last_b_index = var_set.get_index(var_names::b, m, max_K);
    //     const uint v0 = 0, vf = 0;
    //     const Eigen::Vector3d fv0 = traj.getVel(0, t_list[0][0]);
    //     const Eigen::Vector3d fvf = traj.getVel(m, t_list[m][max_K]);
    //     // TODO: setting to something other than zero requires fv0 and fvf
    //     // otherwise, we can just set b0 and bf to zero
    //     for (uint i = 0; i < 3; i++) {

    //     }
    //     A(cur_A_index, init_b_index) = 1;
    // }


    // ------------------------------------------------------------------
    // Debug functions --------------------------------------------------
    // ------------------------------------------------------------------
    // Print matrices
    // save_matrices_to_file (n_variables, size_G, n_eq, n_orthants, n_soc,
    //     soc_size, sG, JG, IG, sA, JA, IA, c, h, b);

    // Print SOCP problem to files (easy to read)
    // save_socp_constraints_to_files (c, A, b, G, h, n_orthants, soc_size, var_set);

    
    // ------------------------------------------------------------------
    // Set ECOS ---------------------------------------------------------
    // // ---------------------------------------------------------------
    pwork *mywork;
    idxint exitflag;
    uint n_exp_cones = 0;
    pfloat *c_ecos = c.data();
    pfloat *h_ecos = h.data();
    pfloat *b_ecos = b.data();
    pfloat *Gx = sG.data(); idxint *Gj = JG.data(); idxint *Gi = IG.data();
    pfloat *Ax = sA.data(); idxint *Aj = JA.data(); idxint *Ai = IA.data();
    idxint *soc_sizes = soc_size.data();
    mywork = ECOS_setup(n_variables, size_G, n_eq, n_orthants, n_soc, 
                        soc_size.data(), n_exp_cones, Gx, Gj, Gi, Ax, Aj, Ai,
                        c_ecos, h_ecos, b_ecos);
    ros::Time t2 = ros::Time::now();
    std::cout << "Setup time: " << (t2 - t1).toSec() << std::endl;

    mywork->stgs->feastol = 1e-6;
    mywork->stgs->abstol = 1e-6;
    mywork->stgs->reltol = 1e-6;
    mywork->stgs->feastol_inacc = 1e-4;
    mywork->stgs->abstol_inacc = 1e-4;
    mywork->stgs->reltol_inacc = 1e-4;
    mywork->stgs->nitref = 9;
    mywork->stgs->gamma = 0.99;
    mywork->stgs->delta = 2e-7;
    mywork->stgs->eps = 1e-9;
    mywork->stgs->maxit = 100;
    mywork->stgs->verbose = 0;


    // std::cout << "setup done" << std::endl;
    if( mywork != NULL ){
        /* solve */
        exitflag = ECOS_solve(mywork); 
    } else {
        exitflag = ECOS_FATAL;
    }

    ros::Time t3 = ros::Time::now();
    std::cout << "Solution time: " << (t3 - t2).toSec() << std::endl;

    if (!check_exit_flag(exitflag)) {
        ECOS_cleanup(mywork, 0);
        return 0;
    }

    // Retrieve 'a' and 'b' solutions
    pfloat* sol_X = mywork->x;

    // Vectors of vectors, splitting segments
    std::vector<std::vector<pfloat>> sol_a;
    sol_a.resize(n_segments);
    idxint idx = 0;
    for (uint i = 0; i <= m; i++) {
        const uint max_Ka  = var_set.a_.segments[i].K_;
        sol_a[i].resize(max_Ka + 1);
        for (uint j = 0; j <= max_Ka; j++) {
            sol_a[i][j] = *(sol_X + idx);
            idx++;
            // std::cout << "(" << i << ", " << j << ") = " << sol_a[i][j] << std::endl;
        }
    }

    std::vector<std::vector<pfloat>> sol_b;
    sol_b.resize(n_segments);
    idx = var_set.b_.initial_index_;
    // std::cout << "sol b: " << std::endl;
    for (uint i = 0; i <= m; i++) {
        const uint max_Kb  = var_set.b_.segments[i].K_;
        sol_b[i].resize(max_Kb + 1);
        for (uint j = 0; j <= max_Kb; j++) {
            sol_b[i][j] = *(sol_X + idx);
            idx++;
            // std::cout << "(" << i << ", " << j << ") = " << sol_b[i][j] << std::endl;
        }
    }

    // std::vector<std::vector<pfloat>> sol_c;
    // sol_c.resize(n_segments);
    // idx = var_set.c_.initial_index_;
    // for (uint i = 0; i <= m; i++) {
    //     const uint max_Kc  = var_set.c_.segments[i].K_;
    //     sol_c[i].resize(max_Kc + 1);
    //     for (uint j = 0; j <= max_Kc; j++) {
    //         sol_c[i][j] = *(sol_X + idx);
    //         idx++;
    //         // std::cout << "(" << i << ", " << j << ") = " << sol_c[i][j] << std::endl;
    //     }
    // }
    // std::cout << "returned values for c" << std::endl;

    // std::vector<std::vector<pfloat>> sol_d;
    // sol_d.resize(n_segments);
    // idx = var_set.d_.initial_index_;
    // for (uint i = 0; i <= m; i++) {
    //     const uint max_Kd  = var_set.d_.segments[i].K_;
    //     sol_d[i].resize(max_Kd + 1);
    //     for (uint j = 0; j <= max_Kd; j++) {
    //         sol_d[i][j] = *(sol_X + idx);
    //         idx++;
    //         // std::cout << "(" << i << ", " << j << ") = " << sol_d[i][j] << std::endl;
    //     }
    // }
    // std::cout << "returned values for d" << std::endl;

    // Check if solution is compliant with constraints
    // check_constraints_compliance (n_variables, n_orthants, n_eq,
    //     sol_X, sol_a, sol_b, sol_c, sol_d, A, b, G, h);


    // Stacking the output results
    int max_K = -1;
    for(int i = 0; i < n_segments; i++) {
        const int cur_K = (int)k_list[i];
        if(cur_K > max_K)
            max_K = cur_K;
    }

    time_allocator = new Allocator(n_segments, new_ds, max_K, maxVel, maxAcc, maxJer_s);
    for(int k = 0; k < n_segments; k++) {
        double T = 0.0;
        int K = (int)k_list[k];

        time_allocator->K(k) = K;
        std::vector<pfloat> a_k   = sol_a[k];
        std::vector<pfloat> b_k   = sol_b[k];

        for(int i = 0; i <= K; i++) {
            if(i <  K) {
                time_allocator->a(k, i) = a_k[i];

                if( b_k[i] <= 0.0 || b_k[i+1] <= 0.0 )
                    T += 0.0;
                else
                    T += 1.0 * 2 * new_ds/(sqrt(b_k[i]) + sqrt(b_k[i+1]));
                
                time_allocator->time(k, i) = T;
                
                if(i == 0)
                    time_allocator->time_acc(k, i) = time_allocator->time(k, i) / 2.0;
                else
                    time_allocator->time_acc(k, i) = (time_allocator->time(k, i) + time_allocator->time(k, i - 1)) / 2.0;
            }    
            
            time_allocator->b(k, i) = b_k[i];
            time_allocator->s(k, i) = t_list[k](i);
        }

    }

    
    /* clean up memory */
    ECOS_cleanup(mywork, 0);


    return 1;
}

bool MinimumTimeOptimizer::check_exit_flag (const idxint &exitflag) {
    switch (exitflag) {
        case ECOS_OPTIMAL:
            ROS_INFO("[ECOS] Optimal solution found!");
            break;
        case ECOS_PINF:
            ROS_WARN("[ECOS] Certificate of primal infeasibility found!");
            return 0;
        case ECOS_DINF:
            ROS_WARN("[ECOS] Certificate of dual infeasibility found!");
            return 0;
        case ECOS_INACC_OFFSET:
            ROS_WARN("[ECOS] Optimal solution found subject to reduced tolerances!");
            break;
        case 11:
            ROS_WARN("[ECOS] Certificate of primal infeasibility found subject to reduced tolerances!");
            return 0;
        case 12:
            ROS_WARN("[ECOS] Certificate of dual infeasibility found subject to reduced tolerances!");
            return 0;
        case ECOS_MAXIT:
            ROS_WARN("[ECOS] Maximum number of iterations reached!");
            return 0;
        case ECOS_NUMERICS:
            ROS_WARN("[ECOS] Numerical problems (unreliable search direction)!");
            return 0;
        case ECOS_OUTCONE:
            ROS_WARN("[ECOS] Numerical problems (slacks or multipliers outside cone)!");
            return 0;
        case ECOS_SIGINT:
            ROS_WARN("[ECOS] Interrupted by signal or CTRL-C!");
            return 0;
        case ECOS_FATAL:
            ROS_WARN("[ECOS] Unknown problem in solver!");
            return 0;
    }
    return 1;
}

void MinimumTimeOptimizer::vector_to_string (
        const std::vector<idxint> &vec, const std::string &type,
        const std::string &var_name, std::string *out) {
    *out = type + " " + var_name + "[" + std::to_string(vec.size()) + "] = {";
    for (uint i = 0; i < vec.size(); i++) {
        out->append(std::to_string(vec[i]));
        if (i != vec.size() - 1) {
            out->append(", ");
        }
    }
    out->append("};\n");
}

void MinimumTimeOptimizer::vector_to_string(
        const std::vector<double> &vec, const std::string &type,
        const std::string &var_name, std::string *out) {
    *out = type + " " + var_name + "[" + std::to_string(vec.size()) + "] = {";
    for (uint i = 0; i < vec.size(); i++) {
        out->append(std::to_string(vec[i]));
        if (i != vec.size() - 1) {
            out->append(", ");
        }
    }
    out->append("};\n");
}

void MinimumTimeOptimizer::vector_to_string(
        const Eigen::VectorXd &vec, const std::string &type,
        const std::string &var_name, std::string *out) {
    *out = type + " " + var_name + "[" + std::to_string(vec.size()) + "] = {";
    for (uint i = 0; i < vec.size(); i++) {
        out->append(std::to_string(vec[i]));
        if (i != vec.size() - 1) {
            out->append(", ");
        }
    }
    out->append("};\n");
}

void MinimumTimeOptimizer::save_matrices_to_file (
    const uint &n_variables, const uint &size_G, const uint &n_eq,
    const uint &n_orthants, const uint &n_soc, const std::vector<idxint> &soc_size,
    const std::vector<double> &sG, const std::vector<idxint> &JG, const std::vector<idxint> &IG,
    const std::vector<double> &sA, const std::vector<idxint> &JA, const std::vector<idxint> &IA,
    const Eigen::VectorXd &c, const Eigen::VectorXd &h, const Eigen::VectorXd &b) {
    ofstream myfile;
    myfile.open ("socp_matrices.txt");
    myfile << "idxint n_var = " << n_variables << ";" << std::endl;
    myfile << "idxint size_G = " << size_G << ";" <<  std::endl;
    myfile << "idxint n_eq = " << n_eq << ";" <<  std::endl;
    myfile << "idxint n_orthants = " << n_orthants << ";" <<  std::endl;
    myfile << "idxint n_soc = " << n_soc << ";"  << ";" <<  std::endl;
    std::string str_vec;

    vector_to_string(soc_size, "idxint", "soc_size", &str_vec);
    myfile << str_vec;
    vector_to_string(sG, "idxint", "Gpr", &str_vec);
    myfile << str_vec;
    vector_to_string(JG, "idxint", "Gjc", &str_vec);
    myfile << str_vec;
    vector_to_string(IG, "idxint", "Gir", &str_vec);
    myfile << str_vec;
    vector_to_string(sA, "idxint", "Apr", &str_vec);
    myfile << str_vec;
    vector_to_string(JA, "idxint", "Ajc", &str_vec);
    myfile << str_vec;
    vector_to_string(IA, "idxint", "Air", &str_vec);
    myfile << str_vec;
    vector_to_string(c, "idxint", "c", &str_vec);
    myfile << str_vec;
    vector_to_string(h, "idxint", "h", &str_vec);
    myfile << str_vec;
    vector_to_string(b, "idxint", "b", &str_vec);
    myfile << str_vec;

    myfile.close();
}

void MinimumTimeOptimizer::save_socp_constraints_to_files (
        const Eigen::VectorXd &c, const Eigen::MatrixXd &A,
        const Eigen::VectorXd &b, const Eigen::MatrixXd &G,
        const Eigen::VectorXd &h, const uint &n_orthants,
        const std::vector<idxint> &soc_size,
        const ecos_sol::variable_set &var_set) {
    uint size_G = G.rows();

    // Print details about the problem
    ofstream costfile, eq_file, ineq_file, soc_file;
    costfile.open ("socp_cost.txt");
    bool first = true;
    std::string var_str;
    costfile << "Cost function:\n";
    for (uint i = 0; i < c.size(); i++) {
        if (first && (c[i] != 0)) {
            var_set.index_to_var_ik(i, &var_str);
            costfile << c[i];
            costfile << "." + var_str;
            first = false;
        } else if (c[i] != 0) {
            var_set.index_to_var_ik(i, &var_str);
            costfile << " + ";
            costfile << c[i];
            costfile << "." + var_str;
        }
    }
    costfile.close();

    eq_file.open ("socp_equalities.txt");
    eq_file << "\n\nEquality constraints:\n";
    first = true;
    for (uint i = 0; i < A.rows(); i++) {
        for (uint j = 0; j < A.cols(); j++) {
            if (first && (A(i,j) != 0)) {
                var_set.index_to_var_ik(j, &var_str);
                eq_file << A(i,j);
                eq_file << "." + var_str;
                first = false;
            } else if (A(i,j) != 0) {
                var_set.index_to_var_ik(j, &var_str);
                if (A(i,j) >= 0) {
                    eq_file << " + ";
                } else {
                    eq_file << " - ";
                }
                eq_file << std::abs(A(i,j));
                eq_file << "." + var_str;
            }
        }
        eq_file << " = ";
        eq_file << b(i);
        eq_file << "\n";
    }
    eq_file.close();
    
    ineq_file.open ("socp_inequalities.txt");
    ineq_file << "\n\nInequality constraints:\n";
    for (uint i = 0; i < n_orthants; i++) {
        first = true;
        for (uint j = 0; j < G.cols(); j++) {
            if (first && (G(i,j) != 0)) {
                var_set.index_to_var_ik(j, &var_str);
                ineq_file << G(i,j);
                ineq_file << "." + var_str;
                first = false;
            } else if (G(i,j) != 0) {
                var_set.index_to_var_ik(j, &var_str);
                if (G(i,j) >= 0) {
                    ineq_file << " + ";
                } else {
                    ineq_file << " - ";
                }
                ineq_file << std::abs(G(i,j));
                ineq_file << "." + var_str;
            }
        }
        ineq_file << " <= ";
        ineq_file << h(i);
        ineq_file << "\n";
    }
    ineq_file.close();

    soc_file.open ("socp_soc_constraints.txt");
    soc_file << "\n\nSOC constraints (skip two lines per cone):\n";
    uint cur_index_soc = 0;
    uint cur_soc = 0;
    for (uint i = n_orthants; i < size_G; i++) {
        first = true;
        for (uint j = 0; j < G.cols(); j++) {
            if (first && (G(i,j) != 0)) {
                var_set.index_to_var_ik(j, &var_str);
                soc_file << G(i,j);
                soc_file << "." + var_str;
                first = false;
            } else if (G(i,j) != 0) {
                var_set.index_to_var_ik(j, &var_str);
                if (G(i,j) >= 0) {
                    soc_file << " + ";
                } else {
                    soc_file << " - ";
                }
                soc_file << std::abs(G(i,j));
                soc_file << "." + var_str;
            }
        }
        if (h(i) >= 0) {
            soc_file << " + ";
        } else {
            soc_file << " - ";
        }
        soc_file << std::abs(h(i));
        soc_file << "\n";

        // Checks if we got to a new SOC
        cur_index_soc++;
        if (soc_size[cur_soc] == cur_index_soc) {
            cur_soc++;
            cur_index_soc = 0;
            soc_file << "\n";
        }
    }
}

void MinimumTimeOptimizer::check_constraints_compliance (
        const uint &n_variables, const uint &n_orthants, 
        const uint &n_eq, const pfloat *sol_X,
        const std::vector<std::vector<pfloat>> &sol_a,
        const std::vector<std::vector<pfloat>> &sol_b,
        const std::vector<std::vector<pfloat>> &sol_c,
        const std::vector<std::vector<pfloat>> &sol_d,
        const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
        const Eigen::MatrixXd &G, const Eigen::VectorXd &h) {
    // Get the whole solution into an Eigen vector
    Eigen::VectorXd solution(n_variables);
    for (uint i = 0; i < n_variables; i++) {
        solution[i] = *(sol_X + i);
    }

    // Check for equalities
    for (uint i = 0; i < n_eq; i++) {
        const double equality = A.row(i)*solution - b(i);
        if (std::fabs(equality) > 0.00000001) {
            std::cout << "Equality not satisfied in line" << 
                         i << ": " << equality << std::endl;
        }
    }
    std::cout << "Checked for equalities!" << std::endl;

    // Check for orthants
    for (uint i = 0; i < n_orthants; i++) {
        const double inequality = h(i) - G.row(i)*solution;
        if (inequality < 0) {
            std::cout << "Inequality not satisfied in line " << 
                          i << ": " << inequality << std::endl;
        }
    }
    std::cout << "Checked for inequalities!" << std::endl;

    // Check if the rotated cone is being satisfied
    for(uint i = 0; i < sol_d.size(); i++) {
        for (uint j = 0; j < sol_d[i].size(); j++) {
            double lhs = sol_d[i][j]*(sol_c[i][j] + sol_c[i][j+1]);
            if (std::fabs(lhs - 1) >= 0.000001) {
                std::cout << "(" << i << ", " << j << ") = " << lhs << std::endl;
            }
        }
    }
    std::cout << "Checked for rotated cone!" << std::endl;

    // Check if the soc is being satisfied
    for(uint i = 0; i < sol_b.size(); i++) {
        for (uint j = 0; j < sol_b[i].size(); j++) {
            double lhs = pow(sol_b[i][j] + 1, 2);
            double rhs = pow(sol_b[i][j] - 1, 2) + pow(2*sol_c[i][j], 2);
            if (std::fabs(lhs - rhs) >= 0.000001) {
                std::cout << "(" << i << ", " << j << ") = " << lhs - rhs << std::endl;
            }
        }
    }
    std::cout << "Checked for second order cone!" << std::endl;
}

}  // namespace ecos_sol