#include <fdaPDE/core/fdaPDE/core.h>
#include <fdaPDE/src/solvers/utility.h>
#include "../include/utils.h"
#include <unsupported/Eigen/SparseExtra>
#include <string>
#include <vector>

using namespace fdapde;

int main (int argc, char* argv[]){
    
    using vector_t = Eigen::Matrix<double, Dynamic, 1>; 
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>; 
    
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using sparsesolver_t = Eigen::SparseLU<sparse_matrix_t>;
    constexpr int local_dim = 2;
    using PointT = Eigen::Matrix<double, local_dim, 1>;
    
    int N = 11;
    int M = 11;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--N" && i + 1 < argc) { 
            N = std::stoi( argv[++i] ) ;
        }

    }

    Triangulation<2, 2> unit_square = Triangulation<2, 2>::UnitSquare(N, cache_cells);

    std::string input_dir = "input/N_" + std::to_string(N) + "/"; 
    std::string mesh_dir = input_dir + "mesh/";
    
    std::string command_str = "mkdir -p " + mesh_dir; 
    system(command_str.c_str());
    Eigen::saveMarket(unit_square.nodes(), mesh_dir + "points.mtx");
    Eigen::saveMarket(unit_square.cells(), mesh_dir + "cells.mtx");
    Eigen::saveMarket(unit_square.boundary_nodes().as_eigen_matrix(), mesh_dir + "boundary.mtx");

    // FE space
    auto quad = P1<1>;
    
    FeSpace Vh(unit_square, quad);   // piecewise linear continuous scalar finite elements
    TrialFunction u(Vh);
    TestFunction  v(Vh);
    
    auto& dof_handler = Vh.dof_handler();
    int n_dofs = dof_handler.n_dofs();
    matrix_t dofs =  unit_square.nodes(); // fe_dofs_coords(unit_square, quad);
    std::cout << dofs.rows() << " " << dofs.cols() << std::endl;
    unit_square.mark_boundary(/* as = */ 0);
    unit_square.mark_boundary(/* as = */ 4, /* where = */ [](const typename Triangulation<2, 2>::EdgeType& e) {
       return (e.node(0)[0] == 0 && e.node(1)[0] == 0);
    });
    //dof_handler.set_dirichlet_constraint(/* on = */ 0, /* data = */ g_D);
    
    FeFunction u_prev(Vh);

    auto coe = [](double x) -> double {
        return 0.5*std::sin(5.*std::numbers::pi*x)*std::exp(-x*x) + 1.;
    };

    double tf = 1.0;
    double t = 0.0;
    
    auto vicini = [&tf, &t, &coe](const PointT& p) {
        double a = (9./(5.*tf))*t - 2.;
        return abs( std::sin( 2.*std::numbers::pi*(  coe(p[1])*p[0]*cos(a) - p[1]*std::sin(a) ) ) * 
                    std::cos( 2.*std::numbers::pi*(  coe(p[1])*p[0]*cos(a + std::numbers::pi/2.) + coe(p[0])*p[1]*std::sin( (std::numbers::pi/2.)*a ) ) )
        );
    };

    ScalarField<2> ic;
    ic = vicini;
     
    double DeltaT = 1e-1;
    int n_times = M+1;
    double T = (n_times-1) * DeltaT;
    vector_t time_mesh = vector_t::Zero(n_times);
    for(int t = 0; t < n_times; t++) time_mesh[t] = DeltaT*t;
    
    std::string data_dir = input_dir;
    Eigen::saveMarket(time_mesh.head(n_times-1), input_dir + "time_mesh.mtx");
    

    vector_t IC = vector_t::Zero(n_dofs);
    for(int i = 0; i < n_dofs; ++i) IC[i] = ic(dofs.row(i));
    std::cout << dofs.rows() << " " << dofs.cols() << std::endl;
    std::cout << unit_square.nodes().rows() << " " << unit_square.nodes().cols() << std::endl;
    
    std::cout << IC.minCoeff() << " " << IC.maxCoeff() << std::endl;
        
    u_prev = IC;

    // laplacian operator bilinear form
    double mu = 1e-2; // 0.01
    std::vector<double> alphas{1.0}; 

    Eigen::Matrix<double, 2,2> K;
    K << 1.1*mu, 0, 0, 0.9*mu;
   
    for (int r=0; r<alphas.size(); ++r){

        double alpha = alphas[r];

        //double alpha = alpha_mtx(k,i);
        auto a = integral(unit_square)(dot(K*grad(u), grad(v)));
        auto reac= integral(unit_square)(u_prev*u*v);
        auto mass = integral(unit_square)(u*v);

        matrix_t solution = matrix_t::Zero(n_dofs, n_times);
        matrix_t diffusion = matrix_t::Zero(n_dofs, n_times);
        
        solution.col(0) = IC;
        diffusion.col(0) = IC;

        auto M = mass.assemble();
        auto A = a.assemble();
        vector_t rhs = vector_t::Zero(n_dofs);
        
        for(int t = 0; t < n_times-1; ++t){
            u_prev = solution.col(t);
            auto R = reac.assemble();
            sparse_matrix_t S = (1./DeltaT*M + A) - alpha*(M - R);
            
            rhs = 1./DeltaT * M*solution.col(t);
            
            sparsesolver_t lin_solver;
            lin_solver.analyzePattern(S); 
            lin_solver.factorize(S);
            solution.col(t+1) = lin_solver.solve(rhs);
            
            //---
            S = (1./DeltaT*M + A);
            rhs = 1./DeltaT * M*diffusion.col(t);
            lin_solver.analyzePattern(S); 
            lin_solver.factorize(S);
            diffusion.col(t+1) = lin_solver.solve(rhs);  
        }

        matrix_t test_locs = expand_grid(std::vector<vector_t>{
                                         vector_t::LinSpaced(60, 0.05, 0.95), 
                                         vector_t::LinSpaced(60, 0.05, 0.95)});
        auto psi_test = internals::point_basis_eval(Vh, test_locs);
    
        matrix_t test_obs = matrix_t::Zero(test_locs.rows(), n_times-1); 
        for(int t = 1; t<n_times; ++t) test_obs.col(t-1) = psi_test * solution.col(t);
        
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(2) << alpha;

        std::string output_dir = data_dir + ss.str() + "/";
        command_str = "mkdir -p " + output_dir; 
        std::cout << "Save data in      " << output_dir << std::endl;
        system(command_str.c_str());
    
        Eigen::saveMarket(solution.rightCols(n_times - 1), output_dir + "fisher_kpp.mtx");
        Eigen::saveMarket(diffusion.rightCols(n_times - 1), output_dir + "diffusion.mtx");
        
        Eigen::saveMarket(test_obs.rightCols(n_times - 1), output_dir + "test_obs.mtx");
        Eigen::saveMarket(test_locs, output_dir + "test_locs.mtx");

        test_obs = matrix_t::Zero(test_locs.rows(), n_times-1); 
        for(int t = 1; t<n_times; ++t) test_obs.col(t-1) = psi_test * diffusion.col(t);

        Eigen::saveMarket(test_obs, output_dir + "test_obs-diffusion.mtx");
    }

    // 30 simulation repetitions
    std::vector<int> n_locs{250};
    int n_sim = 30;
    std::mt19937 gen(12345);
    for (int r=0; r<alphas.size(); ++r){
        
        double alpha = alphas[r];
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(2) << alpha;
        std::string output_dir = data_dir + ss.str() + "/";

        auto solution = read_mtx<double>(output_dir + "fisher_kpp.mtx");
        std::cout << "Read 'field' in   " << output_dir << std::endl;
        for(int i=0; i<n_locs.size(); ++i){
            std::string locs_dir = output_dir + std::to_string(n_locs[i]) + "/";
            
            for(int k = 0; k < n_sim; ++k){
                std::string sim_dir = locs_dir + std::to_string(k) + "/";
    
                matrix_t locs = unit_square.sample(n_locs[i]);
            
                auto psi = internals::point_basis_eval(Vh, locs);
                matrix_t obs = matrix_t::Zero(n_locs[i], M); 
            
                for(int t=0; t < M; ++t) obs.col(t) = psi * solution.col(t);
                
                Eigen::saveMarket(obs, sim_dir + "obs_no_noise.mtx");
                obs += noise(n_locs[i], M, 0.05, gen);

                command_str = "mkdir -p " + sim_dir; 
                system(command_str.c_str());
    
                Eigen::saveMarket(locs, sim_dir + "locs.mtx");
                Eigen::saveMarket(obs, sim_dir + "obs.mtx");
                std::cout << "  Save data in   " << sim_dir << std::endl;
            }
        }
    } 

    command_str = "chown -R 1000:1000 " + input_dir + "/"; 
    system(command_str.c_str());
    return 0;
}
