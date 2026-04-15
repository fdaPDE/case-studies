#include <fdaPDE/fdapde.h>
#include "../include/utils.h"
#include "../include/iterative_method.h"
#include "../include/kFoldCV.h"
#include <unsupported/Eigen/SparseExtra>

using namespace fdapde;


int main(int argc, char *argv[]){
    
    using vector_t = Eigen::Matrix<double, Dynamic, 1>; 
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>; 
    
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using sparsesolver_t = Eigen::SparseLU<sparse_matrix_t>;
    constexpr int local_dim = 2;
    using PointT = Eigen::Matrix<double, local_dim, 1>;
    
    std::string data_dir = "input/";
    std::string mesh_dir = "../simulation_1/input/mesh/";
    std::string sim = std::string(argv[1]); // 0, ..., 29
    std::string sim_dir = data_dir  + "/" + sim + "/";

    std::string output_dir = executable_name(argv[0]) + "/" + sim + "/";
    std::string command_str = "mkdir -p " + output_dir;
    system(command_str.c_str());
    auto test_locs = read_mtx<double>("../simulation_1/input/1.00/test_locs.mtx");
    auto nodes = read_mtx<double>(mesh_dir + "points.mtx");
    auto cells = read_mtx<int>(mesh_dir + "cells.mtx");
    auto boundary = read_mtx<int>(mesh_dir + "boundary.mtx");
    
    std::cout << nodes.rows() << " " << nodes.cols() <<  std::endl;
    std::cout << cells.rows() << " " << cells.cols() <<  std::endl;
    std::cout << boundary.rows() << " " << boundary.cols() <<  std::endl;
    Triangulation<2, 2> unit_square = Triangulation<2, 2>(nodes, cells, boundary);

    FeSpace Vh(unit_square, P1<1>);
    TrialFunction u(Vh);
    TestFunction v(Vh);
    ZeroField<2> f;

    matrix_t time_mesh = read_mtx<double>(mesh_dir + "time_mesh.mtx");
    matrix_t time_locs = time_mesh.block(1, 0, time_mesh.size()-1, 1);
    int n_times = time_mesh.rows();
    Triangulation<1, 1> T(time_mesh);

    Eigen::Matrix<int, Dynamic, Dynamic> incidence_matrix = read_csv<int>("input/incidence_matrix_not_unif.csv").as_matrix();
    matrix_t obs = read_mtx<double>(sim_dir + "obs_not_unif.mtx");
    std::cout << incidence_matrix.rows() << " " << incidence_matrix.cols() << std::endl;
    std::cout << obs.rows() << " " << obs.cols() << std::endl;
    
    matrix_t exact = read_mtx<double>("../simulation_1/input/1.00/fisher_kpp.mtx"); 
    std::cout << exact.rows() << " " << exact.cols() << std::endl;
    auto bm = BinaryMatrix<Dynamic, Dynamic> (incidence_matrix.rows(), incidence_matrix.cols());
    bm = incidence_matrix;

    const auto &[psi, measure_vec] =
        internals::areal_basis_eval(Vh, bm);

    auto Psi_ = psi;
    auto D_ = measure_vec.asDiagonal();

    double mu = 0.01; // 0.001
    
    int n_dofs = nodes.rows();
    
    auto a = integral(unit_square)(mu * dot(grad(u), grad(v))); 
    auto F = integral(unit_square)(f * v);
    auto mass = integral(unit_square)(u*v).assemble();
    
    vector_t IC = exact.col(0);
    // -----------------------------------------------------------------------------------------------------------
    vector_t lambdas = vector_t::Ones(17);
    lambdas << 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1, 1, 5, 1e1, 5e1, 1e2, 5e2, 1e3;
    
    vector_t reactions = vector_t::Ones(9);
    reactions << 0., 0.1, 0.25, 0.5, 0.9, 1.0, 1.1, 1.5, 2.0;

    // ---------------------------------------------------------------------------------------------------------------------------
     Eigen::Matrix<double, Dynamic,2> grid = expand_grid( std::vector<vector_t>{lambdas, reactions} );
    Eigen::Matrix<double, Dynamic,2> grid_linear = expand_grid( std::vector<vector_t>{lambdas, vector_t::Zero(1)} );
    
    std::cout << grid.transpose() << std::endl;
    std::cout << " linear \n" << grid_linear.transpose() << std::endl;
    
    Eigen::saveMarket(grid, output_dir + "cv_grids.mtx");
    int K = 10; 

    KFoldCV cv_(incidence_matrix, obs, K);
    
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

            auto incidence_mtx = cv_.X_train(k);
         
            matrix_t obs = cv_.Y_train(k);
             
            auto bm = BinaryMatrix<Dynamic, Dynamic> (incidence_mtx.rows(), incidence_mtx.cols());
            bm = incidence_mtx;

            vector_t  result = vector_t::Zero((n_times-1)*n_dofs);
            if( alpha > 0){

                vector_t response = obs.rightCols(n_times - 1).reshaped();
                GeoFrame data(unit_square, T); 
                auto &l = data.insert_scalar_layer<POLYGON, POINT>(
                "layer", std::pair{bm, time_locs}); 
                l.load_vec("y", response.reshaped() );
                
                SRPDE parabolic("y ~ f", data, fe_ls_parabolic_mono(std::pair{a, F}, IC));
                parabolic.fit(std::vector<double>{lambda, 1.0});
                
                auto lbfgs = LBFGS<Dynamic>(100, 1e-7, 1e-2, 10);
                auto model = iterative_method(unit_square, mu, alpha, bm, time_mesh,
                                 obs.reshaped());

                vector_t g_init = vector_t::Zero(n_dofs*n_times);
                g_init.head(n_dofs * (n_times-1)) = parabolic.misfit();
                model.set_state_initial_condition(IC);
                model.set_g_initial_guess(g_init);
                model.solve(lambda, lbfgs); // BacktrackingLineSearch()
            
                result = model.y().tail((n_times-1)*n_dofs);
                std::cout << "#iters " << model.n_iter() << std::endl;
            }else{
                vector_t response = obs.rightCols(n_times - 1).reshaped();
                GeoFrame data(unit_square, T); 
                auto &l = data.insert_scalar_layer<POLYGON, POINT>(
                "layer", std::pair{bm, time_locs}); 
                l.load_vec("y", response.reshaped() );
                
                SRPDE parabolic("y ~ f", data, fe_ls_parabolic_mono(std::pair{a, F}, IC));
                parabolic.fit(std::vector<double>{lambda, 1.0});
                result = parabolic.f();
            }

            auto test_locs = cv_.X_test(k);
            matrix_t test_obs = cv_.Y_test(k).rightCols(n_times-1);
            auto test_bm = BinaryMatrix<Dynamic, Dynamic> (test_locs.rows(), test_locs.cols()); // 1 x ...
            test_bm = test_locs; 
            auto [psi_test, measure_vect] = internals::areal_basis_eval(Vh, test_bm);

            matrix_t test_vals = matrix_t::Zero(test_locs.rows(), n_times - 1); // t0 butto via
            
            sparse_matrix_t Psi_test = kronecker(Im,psi_test);
            test_vals = Psi_test * result;
            
            res += std::sqrt( (test_vals - test_obs.reshaped()).array().square().mean());
        }
        
        std::cout << "lambda " << lambda << " alpha " << alpha <<  " cv_error = " << 1./cv_.k() * res<<  std::endl;
        return 1./cv_.k() * res;

    };

    GridSearch<2> optimizer_linear;
    optimizer_linear.optimize(SSE, grid_linear); 
    
    GridSearch<2> optimizer;
    optimizer.optimize(SSE, grid); 
    
    vector_t rmse = vector_t::Ones(1);
    matrix_t test_obs = read_mtx<double>("../simulation_1/input/1.00/test_obs.mtx");  
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

    std::cout << "opt (linear) " << optimizer.optimum()[0] << std::endl;
    std::cout << "values (linear) " << values_linear << std::endl;
    // ---- output kFold 

    vector_t response = obs.rightCols(n_times - 1).reshaped();
    GeoFrame data(unit_square, T); 
    auto &l = data.insert_scalar_layer<POLYGON, POINT>(
                "layer", std::pair{bm, time_locs}); 
    l.load_vec("y", response.reshaped() );   
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
    auto model = iterative_method(unit_square, mu, optimizer.optimum()[1], bm, time_mesh, obs.reshaped());
    vector_t g_init = vector_t::Zero(n_dofs*n_times);
    g_init.head(n_dofs * (n_times -1 )) = parabolic.misfit();
    model.set_g_initial_guess(g_init);
    model.set_state_initial_condition(IC);
    model.solve(optimizer.optimum()[0]); 
    vector_t result = model.y().tail((n_times-1)*n_dofs);
    
    test_vals = Psi_test * result; 
    rmse[0] = std::sqrt( (test_vals - test_obs.reshaped()).array().square().mean());

    tmp.tail((n_times-1)*n_dofs) = result;
    std::cout << "rmse (NL)" << rmse << std::endl;
    Eigen::saveMarket(rmse, output_dir + "rmse_iterative.mtx");
    Eigen::saveMarket(tmp, output_dir + "estimate_iterative.mtx");
    Eigen::saveMarket(model.u(), output_dir + "misfit_iterative.mtx");
    
    command_str = "chown -R 1000:1000 " + output_dir;
    system(command_str.c_str());
}
