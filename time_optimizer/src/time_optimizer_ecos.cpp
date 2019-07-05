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
    std::cout << "segment times: " << _T.transpose() << std::endl;

    // int num_a[n_segments], num_b[n_segments], num_c[n_segments], num_d[n_segments];

    std::vector<Eigen::VectorXd> t_list;
    std::vector<uint> k_list;
    uint num_X = 0;
    uint sum_K = 0;
    double new_ds; // d_s might change slightly if duration/d_s is not an integer
    for(int j = 0; j < n_segments; j++) {
        double duration  = _T(j);
        uint K = ceil(duration / d_s);
        // std::cout << K << std::endl;

        // Get lists of times
        sum_K = sum_K + K;
        k_list.push_back(K);
        Eigen::VectorXd t_k(K + 1);
        new_ds = duration/K;
        for(int i = 0; i < K + 1; i++) {
            t_k(i) = i * new_ds;
        }
        t_list.push_back(t_k);

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
    uint n_eq   = n_vel_continuity + n_b2a;
    std::cout << "Number of inequalities: " << n_ineq << std::endl;
    std::cout << "Number of equalities: "   <<   n_eq << std::endl;
    std::cout << "Number of variables: "    << n_variables << std::endl;

    Eigen::VectorXd c(n_variables);
    Eigen::MatrixXd G(size_G, n_variables);
    Eigen::VectorXd h(size_G);
    Eigen::MatrixXd A(n_eq, n_variables);
    Eigen::VectorXd b(n_eq);

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
            G(cur_G_index, bk_index) = -1;
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
        // std::cout << "maxk b: " << max_K << std::endl;
        Eigen::VectorXd time = t_list[i];
        for (uint k = 0; k <= max_K; k++) {
            const uint bk_index = var_set.get_index(var_names::b, i, k);
            const Eigen::Vector3d f_prime = traj.getVel(i, time(k));

            // Add for x, y, and z directions
            for (uint idx = 0; idx < 3; idx++) {
                G(cur_G_index, bk_index) = f_prime(idx)*f_prime(idx);
                h(cur_G_index) = max_vel_sqr;
                cur_G_index++;
                n_orthants++;
            }

            // G(cur_G_index, bk_index) = -1;
            // h(cur_G_index) = 0;
            // std::cout << var_set.get_index(var_names::b, i, j) << std::endl;
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
                G(cur_G_index, ak_index) =          f_prime(idx);
                G(cur_G_index, bk_index) =      0.5*f_prime2(idx);
                G(cur_G_index, bkplus1_index) = 0.5*f_prime2(idx);
                h(cur_G_index) = maxAcc;
                cur_G_index++;
                n_orthants++;
            }

            // Add for x, y, and z directions (lower bound)
            for (uint idx = 0; idx < 3; idx++) {
                G(cur_G_index, ak_index) =          -f_prime(idx);
                G(cur_G_index, bk_index) =      -0.5*f_prime2(idx);
                G(cur_G_index, bkplus1_index) = -0.5*f_prime2(idx);
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
                G(cur_G_index, ak_index) =           f_prime(idx);
                G(cur_G_index, akplus1_index) =     -f_plus1_prime(idx);
                G(cur_G_index, bk_index) =       0.5*f_prime2(idx);
                G(cur_G_index, bkplus1_index) =  0.5*f_prime2(idx) - 0.5*f_plus1_prime2(idx);
                G(cur_G_index, bkplus2_index) = -0.5*f_plus1_prime2(idx);
                h(cur_G_index) = maxdAcc;
                cur_G_index++;
                n_orthants++;
            }

            // Add for x, y, and z directions (lower bound)
            for (uint idx = 0; idx < 3; idx++) {
                G(cur_G_index, ak_index) =          -f_prime(idx);
                G(cur_G_index, akplus1_index) =      f_plus1_prime(idx);
                G(cur_G_index, bk_index) =      -0.5*f_prime2(idx);
                G(cur_G_index, bkplus1_index) = -0.5*f_prime2(idx) + 0.5*f_plus1_prime2(idx);
                G(cur_G_index, bkplus2_index) =  0.5*f_plus1_prime2(idx);
                h(cur_G_index) = maxdAcc;
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
        const Eigen::Vector3d f_minus1_prime = traj.getVel(i, eval_time2);
        const Eigen::Vector3d f_minus1_prime2 = traj.getAcc(i, eval_time2);


        
        const uint a_iminus1_kminus1_index = var_set.get_index(var_names::a, i-1, max_Ka);
        const uint a_i_0_index = var_set.get_index(var_names::a, i, 0);
        const uint b_iminus1_kminus1_index = var_set.get_index(var_names::b, i-1, max_Ka);
        const uint b_iminus1_k_index = b_iminus1_kminus1_index + 1;
        const uint b_i_0_index = var_set.get_index(var_names::b, i, 0);
        const uint b_i_1_index = var_set.get_index(var_names::b, i, 1);

        // Add for x, y, and z directions (upper bound)
        for (uint idx = 0; idx < 3; idx++) {
            G(cur_G_index, a_iminus1_kminus1_index) =-f_minus1_prime(idx);
            G(cur_G_index, a_i_0_index) =             f_prime(idx);
            G(cur_G_index, b_iminus1_kminus1_index) =-0.5*f_minus1_prime2(idx);
            G(cur_G_index, b_iminus1_k_index) =      -0.5*f_minus1_prime2(idx);
            G(cur_G_index, b_i_0_index) =             0.5*f_prime2(idx);
            G(cur_G_index, b_i_1_index) =             0.5*f_prime2(idx);
            h(cur_G_index) = maxdAcc;
            cur_G_index++;
            n_orthants++;
        }

        // Add for x, y, and z directions (lower bound)
        for (uint idx = 0; idx < 3; idx++) {
            G(cur_G_index, a_iminus1_kminus1_index) = f_minus1_prime(idx);
            G(cur_G_index, a_i_0_index) =            -f_prime(idx);
            G(cur_G_index, b_iminus1_kminus1_index) = 0.5*f_minus1_prime2(idx);
            G(cur_G_index, b_iminus1_k_index) =       0.5*f_minus1_prime2(idx);
            G(cur_G_index, b_i_0_index) =            -0.5*f_prime2(idx);
            G(cur_G_index, b_i_1_index) =            -0.5*f_prime2(idx);
            h(cur_G_index) = maxdAcc;
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
    const Eigen::Vector3d f_prime_final  = traj.getVel(0, final_time - new_ds);
    const Eigen::Vector3d f_prime2_final = traj.getAcc(0, final_time - new_ds);
    const double init_acc = 0, final_acc = 0;

    // Initial boundary conditions for acceleration (upper bound)
    for (uint idx = 0; idx < 3; idx++) {
        G(cur_G_index, a0_index) = f_prime_init(idx);
        G(cur_G_index, b0_index) = 0.5*f_prime2_init(idx);
        G(cur_G_index, b1_index) = 0.5*f_prime2_init(idx);
        h(cur_G_index) = maxdAcc + init_acc;
        cur_G_index++;
        n_orthants++;
    }

    // Initial boundary conditions for acceleration (lower bound)
    for (uint idx = 0; idx < 3; idx++) {
        G(cur_G_index, a0_index) = -f_prime_init(idx);
        G(cur_G_index, b0_index) = -0.5*f_prime2_init(idx);
        G(cur_G_index, b1_index) = -0.5*f_prime2_init(idx);
        h(cur_G_index) = maxdAcc - init_acc;
        cur_G_index++;
        n_orthants++;
    }

    // Final boundary conditions for acceleration (upper bound)
    for (uint idx = 0; idx < 3; idx++) {
        G(cur_G_index, last_a_index) =   f_prime_final(idx);
        G(cur_G_index, b4last_b_index) = 0.5*f_prime2_final(idx);
        G(cur_G_index, last_b_index) =   0.5*f_prime2_final(idx);
        h(cur_G_index) = maxdAcc + init_acc;
        cur_G_index++;
        n_orthants++;
    }

    // Final boundary conditions for acceleration (lower bound)
    for (uint idx = 0; idx < 3; idx++) {
        G(cur_G_index, last_a_index) =   -f_prime_final(idx);
        G(cur_G_index, b4last_b_index) = -0.5*f_prime2_final(idx);
        G(cur_G_index, last_b_index) =   -0.5*f_prime2_final(idx);
        h(cur_G_index) = maxdAcc - init_acc;
        cur_G_index++;
        n_orthants++;
    }

    // ------------------------------------------------------------------
    // Set second order cones--------------------------------------------
    // ------------------------------------------------------------------
    std::vector<uint> soc_size;

    // Setup the rotated cone constraints -----------------------------
    // dk*(c_{k+1} + c_k)
    for (uint i = 0; i <= m; i++) {
        uint max_Kd = var_set.d_.segments[i].K_;
        for (uint k = 0; k <= max_Kd; k++) {
            const uint dk_index = var_set.get_index(var_names::d, i, k);
            const uint ck_index = var_set.get_index(var_names::c, i, k);
            const uint ckplus1_index = ck_index + 1;

            G(cur_G_index + 1, ck_index) =      1;
            G(cur_G_index + 1, ckplus1_index) = 1;
            G(cur_G_index + 1, dk_index) =     -1;
            G(cur_G_index + 2, ck_index) =     -1;
            G(cur_G_index + 2, ckplus1_index) =-1;
            G(cur_G_index + 2, dk_index) =     -1;
            h(cur_G_index + 0) = 2;
            h(cur_G_index + 1) = 0;
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
            const uint bk_index = var_set.get_index(var_names::b, i, k);
            const uint ck_index = var_set.get_index(var_names::c, i, k);

            G(cur_G_index + 0, bk_index) = -1;
            G(cur_G_index + 0, ck_index) =  0;
            G(cur_G_index + 1, bk_index) =  0;
            G(cur_G_index + 1, ck_index) = -2;
            G(cur_G_index + 2, bk_index) = -1;
            G(cur_G_index + 2, ck_index) =  0;
            h(cur_G_index + 0) = -1;
            h(cur_G_index + 1) =  0;
            h(cur_G_index + 2) =  1;
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
        const uint max_K = var_set.b_.segments[i].K_;
        const uint b_i_K_index = var_set.get_index(var_names::b, i, max_K);
        const uint b_iplus1_0_index = var_set.get_index(var_names::b, i+1, 0);
        A(cur_A_index, b_i_K_index) = 1;
        A(cur_A_index, b_iplus1_0_index) = -1;
        b(cur_A_index) = 0;
        cur_A_index++;
    }

    // Relation between b and a
    // bi_k1 - bi_k - 2.ds.ai_k = 0
    for (uint i = 0; i <= m; i++) {
        uint max_Ka = var_set.a_.segments[i].K_;
        for (uint k = 0; k <= max_Ka; k++) {
            const uint ak_index = var_set.get_index(var_names::a, i, k);
            const uint bk_index = var_set.get_index(var_names::b, i, k);
            const uint bk1_index = bk_index + 1;
            A(cur_A_index, ak_index) = -2*new_ds;
            A(cur_A_index, bk_index) = -1;
            A(cur_A_index, bk1_index) = 1;
            b(cur_A_index) = 0;
            cur_A_index++;
        }
    }

    std::cout << "added inequalities: " << n_soc+n_orthants << std::endl;
    std::cout << "added equalities: "   << cur_A_index      << std::endl;

    ros::Time t0 = ros::Time::now();
    std::vector<double> sG, sA;
    std::vector<int> IG, JG, IA, JA;
    make_sparse_ccs(G, &sG, &IG, &JG);
    make_sparse_ccs(G, &sA, &IA, &JA);
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
    // Set ECOS ---------------------------------------------------------
    // ------------------------------------------------------------------
    


    return 1;
}

}  // namespace scs_sol