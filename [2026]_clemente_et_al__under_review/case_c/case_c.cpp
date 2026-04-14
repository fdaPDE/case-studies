// osservazioni spaziali fisse (sarebbe pronto per gestire anche diversi dati)
// parametro di reazione variabile r = 1.00, n = 125, 250, ..., 2000

#include <fdaPDE/fdapde.h>
#include <fdaPDE/src/solvers/utility.h>
#include "../include/utils.h"
#include "../include/iterative_method.h"
#include "../include/kFoldCV.h"
#include <unsupported/Eigen/SparseExtra>
#include <string>
#include <vector>

using namespace fdapde;
int main(int argc, char *argv[]){
    
    using vector_t = Eigen::Matrix<double, Dynamic, 1>; 
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>; 
    
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using sparsesolver_t = Eigen::SparseLU<sparse_matrix_t>;
    constexpr int local_dim = 2;
    using PointT = Eigen::Matrix<double, local_dim, 1>;
    
   
    std::string mesh_dir = "input/mesh/";
    
    std::string n_locs = std::string( argv[1] );
    std::string sim = std::string(argv[2]); // 0, ..., 29
    std::string r_param = std::string(argv[3]); // 0.20, 0.40, 0.60, 0.80 ...
    std::string M = std::string(argv[4]); // 5, 11, 21, 41 

    std::string data_dir = "input/M_" + M + "/" + r_param + "/";
    std::string sim_dir = data_dir + n_locs  + "/" + sim + "/";


    std::string output_dir = executable_name(argv[0]) + "/M_" +  M + "/" + r_param + "/" + n_locs + "/" + sim + "/";
    std::string command_str = "mkdir -p " + output_dir;
    system(command_str.c_str());
    auto nodes = read_mtx<double>(mesh_dir + "points.mtx");
    auto cells = read_mtx<int>(mesh_dir + "cells.mtx");
    auto boundary = read_mtx<int>(mesh_dir + "boundary.mtx");
    matrix_t test_locs = read_mtx<double>(data_dir + "test_locs.mtx");
    
    Triangulation<2, 2> unit_square = Triangulation<2, 2>(nodes, cells, boundary);

    FeSpace Vh(unit_square, P1<1>);
    TrialFunction u(Vh);
    TestFunction v(Vh);
    ZeroField<2> f;

    matrix_t time_mesh = read_mtx<double>("input/M_" + M + "/time_mesh.mtx");
    matrix_t time_locs = time_mesh.block(1, 0, time_mesh.size()-1, 1);
    int n_times = time_mesh.rows();
    Triangulation<1, 1> T(time_mesh);

    matrix_t locs = read_mtx<double>(sim_dir + "locs.mtx");
    matrix_t obs = read_mtx<double>(sim_dir + "obs.mtx");
    
    matrix_t exact = read_mtx<double>(data_dir + "fisher_kpp.mtx"); 
    //vector_t IC = exact.col(0);

    double mu = 0.01; // 0.001
    
    int n_dofs = nodes.rows();
    //vector_t g_init = vector_t::Random(n_dofs*time_mesh.size()); // + vector_t::Ones(n_dofs*time_locs.size()));
    auto a = integral(unit_square)(mu * dot(grad(u), grad(v))); 
    auto F = integral(unit_square)(f * v);
    auto mass = integral(unit_square)(u*v).assemble();

    vector_t IC = exact.col(0);
    
    vector_t lambdas = vector_t::Ones(9);
    lambdas << 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1, 1;
    
    // vector_t lambdas = vector_t::Ones(1);
    // lambdas << 10;

    vector_t reactions = vector_t::Zero(6); 
    reactions.tail(5) = vector_t::LinSpaced(5, std::stod(r_param) - 0.1, std::stod(r_param) + 0.1);
    
    // vector_t reactions = vector_t::Zero(2); 
    // reactions<< 0. , 1.;
    
    // ---------------------------------------------------------------------------------------------------------------------------
    Eigen::Matrix<double, Dynamic,2> grid = expand_grid( std::vector<vector_t>{lambdas, reactions} );
    Eigen::Matrix<double, Dynamic,2> grid_linear = expand_grid( std::vector<vector_t>{lambdas, vector_t::Zero(1)} );
    
    std::cout << grid.transpose() << std::endl;
    std::cout << " linear \n" << grid_linear.transpose() << std::endl;
    
    Eigen::saveMarket(grid, sim_dir + "cv_grids.mtx");

    int K = 10;
    KFoldCV cv_(locs, obs, K);
    vector_t cv_error = vector_t::Ones(K);

    sparse_matrix_t Im(n_times-1, n_times-1);
    Im.setIdentity();
    sparse_matrix_t Mt = kronecker(Im,mass);
    

  
    ScalarField<2> SSE;
   
    SSE = [&unit_square, &cv_, &time_mesh, &IC, &Vh, &mu, &n_times, &n_dofs, &Im, &T, &a, &F, &time_locs](const vector_t& grid) -> double {
        double lambda = grid[0];
        double alpha = grid[1];
        if (lambda <= 0.0) return std::numeric_limits<double>::max();
        
        double res = 0;
        for(int k=0; k < cv_.k(); ++k){

            matrix_t locs = cv_.X_train(k);
            matrix_t obs = cv_.Y_train(k);
            vector_t  result = vector_t::Zero((n_times-1)*n_dofs);
            if( alpha > 0){

                vector_t response = obs.rightCols(n_times - 1).reshaped();
                GeoFrame data(unit_square, T);
                auto &l = data.insert_scalar_layer<POINT, POINT>("layer", std::pair{locs, time_locs});
                l.load_vec("y", response);
                SRPDE parabolic("y ~ f", data, fe_ls_parabolic_mono(std::pair{a, F}, IC));
                parabolic.fit(std::vector<double>{lambda, 1.0});

                auto model = iterative_method(unit_square, mu, alpha, locs, time_mesh,
                                            obs.reshaped());

                vector_t g_init = vector_t::Zero(n_dofs*n_times);
                g_init.head(n_dofs * (n_times -1 )) = parabolic.misfit();
                model.set_g_initial_guess(g_init);
                model.set_state_initial_condition(IC);
                
                model.solve(lambda); // BacktrackingLineSearch()
                std::cout << "iters: " << model.n_iter() << std::endl;
                matrix_t test_locs = cv_.X_test(k);
                matrix_t test_obs = cv_.Y_test(k).rightCols(n_times-1);
                sparse_matrix_t psi_test = internals::point_basis_eval(Vh, test_locs);
             
                sparse_matrix_t Psi_test = kronecker(Im,psi_test);
            
                result = model.y().tail((n_times-1)*n_dofs);
            
            }else{

                vector_t response = obs.rightCols(n_times - 1).reshaped();
                GeoFrame data(unit_square, T);
                auto &l = data.insert_scalar_layer<POINT, POINT>("layer", std::pair{locs, time_locs});
                l.load_vec("y", response);
                SRPDE parabolic("y ~ f", data, fe_ls_parabolic_mono(std::pair{a, F}, IC));
                parabolic.fit(std::vector<double>{lambda, 1.0});
                result = parabolic.f();
            }

            matrix_t test_locs = cv_.X_test(k);
            matrix_t test_obs = cv_.Y_test(k).rightCols(n_times-1);
            sparse_matrix_t psi_test = internals::point_basis_eval(Vh, test_locs);
            
            sparse_matrix_t Psi_test = kronecker(Im,psi_test);
            vector_t test_vals = Psi_test * result; //;matrix_t::Zero(test_locs.rows(), n_times - 1); // t0 butto via
            res += std::sqrt( (test_vals - test_obs.reshaped()).array().square().mean());
        }
        
        std::cout << "lambda " << lambda << " alpha " << alpha <<  " cv_error = " << 1./cv_.k() * res<<  std::endl;
        return 1./cv_.k() * res;

    };

   
    GridSearch<2> optimizer_linear;
    optimizer_linear.optimize(SSE, grid_linear); 
    
    GridSearch<2> optimizer;
    optimizer.optimize(SSE, grid); 
    
    //
    
    //
    vector_t rmse = vector_t::Ones(1);
    matrix_t test_obs = read_mtx<double>(data_dir + "test_obs.mtx");  
    std::cout << "test_obs " << test_obs.rows() << " " << test_obs.cols() << std::endl;
    auto psi_test = internals::point_basis_eval(Vh, test_locs);
    sparse_matrix_t Psi_test = kronecker(Im, psi_test);

    vector_t values = vector_t::Zero(optimizer.values().size());
    vector_t values_linear = vector_t::Zero(optimizer_linear.values().size());
    
    for(int i = 0; i < optimizer.values().size(); ++i) values[i] = optimizer.values()[i];
    for(int i = 0; i < optimizer_linear.values().size(); ++i) values_linear[i] = optimizer_linear.values()[i];
    
    Eigen::saveMarket(values, output_dir + "cv_errors.mtx");
    Eigen::saveMarket(values_linear, output_dir + "cv_errors_linear.mtx");

    Eigen::saveMarket(optimizer.optimum(), output_dir + "cv_optim.mtx");
    Eigen::saveMarket(optimizer_linear.optimum(), output_dir + "cv_optim_linear.mtx");

    // PARA
    //vector_t values_para = values.head(lambdas.size());
    // Eigen::Index minRow;
    // values_linear.minCoeff(&minRow);
    // double lambda_para = lambdas[minRow]; 
    std::cout << "opt (linear) " << optimizer.optimum()[0] << std::endl;
    std::cout << "values (linear) " << values_linear << std::endl;
    // ---- output kFold 

    vector_t response = obs.rightCols(n_times - 1).reshaped();
    GeoFrame data(unit_square, T);
    auto &l = data.insert_scalar_layer<POINT, POINT>("layer", std::pair{locs, time_locs});
    l.load_vec("y", response);
    SRPDE parabolic("y ~ f", data, fe_ls_parabolic_mono(std::pair{a, F}, IC));

    parabolic.fit( std::vector<double>{ optimizer_linear.optimum()[0], 1.0} );
    vector_t test_vals = Psi_test * parabolic.f();
    rmse[0] = std::sqrt( (test_vals - test_obs.reshaped()).array().square().mean());

    std::cout << "rmse (linear)" << rmse << std::endl;
    Eigen::saveMarket(rmse, output_dir + "rmse_diff_kfold.mtx");
    
    vector_t tmp = vector_t::Zero(n_dofs*n_times);
    tmp.head(n_dofs) = IC;
    tmp.tail(n_dofs*(n_times-1)) = parabolic.f();

    Eigen::saveMarket(tmp, output_dir + "estimate_diff_kfold.mtx");
    Eigen::saveMarket(parabolic.misfit(), output_dir + "misfit_diff_kfold.mtx");

    std::cout << "opt " << optimizer.optimum() << std::endl;    
    std::cout << "values " << optimizer.values() << std::endl;
    
    // nonlinear
    parabolic.fit(std::vector<double>{optimizer.optimum()[0], 1.0});
    auto model = iterative_method(unit_square, mu, optimizer.optimum()[1], locs, time_mesh, obs.reshaped());
    vector_t g_init = vector_t::Zero(n_dofs*n_times);
    g_init.head(n_dofs * (n_times -1 )) = parabolic.misfit();
    model.set_g_initial_guess(g_init);
    model.set_state_initial_condition(IC);
    model.solve(optimizer.optimum()[0]); // BacktrackingLineSearch()
    vector_t result = model.y().tail((n_times-1)*n_dofs);
    
    test_vals = Psi_test * result; //;matrix_t::Zero(test_locs.rows(), n_times - 1); // t0 butto via
    rmse[0] = std::sqrt( (test_vals - test_obs.reshaped()).array().square().mean());

    tmp.tail((n_times-1)*n_dofs) = result;
    std::cout << "rmse (NL)" << rmse << std::endl;
    Eigen::saveMarket(rmse, output_dir + "rmse_iterative.mtx");
    Eigen::saveMarket(tmp, output_dir + "estimate_iterative.mtx");
    Eigen::saveMarket(model.u(), output_dir + "misfit_iterative.mtx");
    
    command_str = "chown -R 1000:1000 " + output_dir;
    system(command_str.c_str());

    return 0;
}