

#include <fdaPDE/core/fdaPDE/core.h>
#include <fdaPDE/src/solvers/utility.h>
#include "../include/utils.h"
#include <unsupported/Eigen/SparseExtra>
#include <string>
#include <vector>

using namespace fdapde;
int main (){
    
    using vector_t = Eigen::Matrix<double, Dynamic, 1>; 
    using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>; 
    
    using sparse_matrix_t = Eigen::SparseMatrix<double>;
    using sparsesolver_t = Eigen::SparseLU<sparse_matrix_t>;
    constexpr int local_dim = 2;
    using PointT = Eigen::Matrix<double, local_dim, 1>;
    
    std::string mesh_dir = "input/mesh/";
    
    auto nodes = read_mtx<double>(mesh_dir + "points.mtx");
    auto cells = read_mtx<int>(mesh_dir + "cells.mtx");
    auto boundary = read_mtx<int>(mesh_dir + "boundary.mtx");
    
    Triangulation<2, 2> unit_square = Triangulation<2, 2>(nodes, cells, boundary);
    FeSpace Vh(unit_square, P1<1>);   // piecewise linear continuous scalar finite elements
    
    std::string command_str;
    vector_t time_mesh = read_mtx<double>(mesh_dir + "time_mesh.mtx");
    int n_times = time_mesh.size();
    std::cout << "time mesh: " << time_mesh.transpose() << std::endl;
    std::cout << "deltaT: " << time_mesh[1] - time_mesh[0] << std::endl; 
    
    std::vector<int> n_locs{125, 250, 500, 1000, 2000}; //{100, 250, 500, 1000, 5000, 10000};

    std::mt19937 gen(12345);
    
    std::vector<double> alphas{0.2, 0.4, 0.6, 0.8, 1.0}; //alphas{1.0};
    //double alpha = 1.;
    int n_sim = 30;

    for (int r=0; r<alphas.size(); ++r){
        
        double alpha = alphas[r];
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(2) << alpha;

        std::string data_dir = "input/" + ss.str() + "/";
        auto solution = read_mtx<double>(data_dir + "fisher_kpp.mtx");
        std::cout << "Read 'field' in   " << data_dir << " rows " << solution.rows() << " cols " << solution.cols() << std::endl;
        for(int i=0; i<n_locs.size(); ++i){
            std::string locs_dir = data_dir + std::to_string(n_locs[i]) + "/";
            for(int k = 0; k < n_sim; ++k){
                std::string sim_dir = locs_dir + std::to_string(k) + "/";
                matrix_t locs = unit_square.sample(n_locs[i]);
            
                auto psi = internals::point_basis_eval(Vh, locs);
                matrix_t obs = matrix_t::Zero(n_locs[i], n_times); 
            
                for(int t=0; t < n_times; ++t) obs.col(t) = psi * solution.col(t);
                
                command_str = "mkdir -p " + sim_dir; 
                system(command_str.c_str());
    
                Eigen::saveMarket(obs, sim_dir + "obs_no_noise.mtx");
                //std::cout << "sd: " << 0.05*std::abs(solution.maxCoeff()) << std::endl;
                obs += noise(n_locs[i], n_times, 0.05*std::abs(solution.maxCoeff()), gen); // minCoeff è zero...

                Eigen::saveMarket(locs, sim_dir + "locs.mtx");
                Eigen::saveMarket(obs, sim_dir + "obs.mtx");
                std::cout << "  Save data in   " << sim_dir << std::endl;
            }
        }
    } 
    command_str = "chown -R 1000:1000 input/"; 
    system(command_str.c_str());

    return 0;
}
