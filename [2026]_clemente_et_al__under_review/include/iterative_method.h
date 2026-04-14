#ifndef __FISHER_KPP_ITERATIVE_H__
#define __FISHER_KPP_ITERATIVE_H__

namespace fdapde {

template<int N>
struct iterative_method {

private:
  void init() {
    FeSpace Vh(domain_, P1<1>);
    TrialFunction uh(Vh);
    TestFunction vh(Vh);

    auto a =
        integral(domain_)(mu_ * dot(grad(uh), grad(vh)) - alpha_ * uh * vh);
    auto mass = integral(domain_)(uh * vh);

    A_ = a.assemble();
    A_.makeCompressed();

    M_ = mass.assemble();
    M_.makeCompressed();

    n_dofs_ = Vh.n_dofs();

    vector_t measure_vec = vector_t::Ones(n_locs_);
    sparse_matrix_t psi(n_locs_,n_dofs_);
    if( incidence_matrix_.size() != 0){
        auto tmp = internals::areal_basis_eval(a.trial_space(), incidence_matrix_);
        Psi_ = std::get<0>(tmp);
        measure_vec = std::get<1>(tmp);

    }else{
        Psi_ = internals::point_basis_eval(a.trial_space(), locations_);
    }

    Psi_.makeCompressed();
    D_ = sparse_matrix_t(n_locs_,n_locs_);
    for(int i = 0; i < n_locs_; ++i) D_.insert(i,i) =  measure_vec[i];

    PsiTDPsi_ = Psi_.transpose() * D_ * Psi_;
    PsiTDPsi_.makeCompressed();


    y_ = vector_t::Zero(n_dofs_ * n_time_locs_);
    u_ = y_;
    p_ = y_;

    sparse_matrix_t I_(n_time_locs_, n_time_locs_);
    I_.setIdentity();

    Mt_ = kronecker(I_, M_);
    Mt_.makeCompressed();

    Psi_t_ = kronecker(I_, Psi_);
    
  }
  
public:
  using vector_t = Eigen::Matrix<double, Dynamic, 1>;
  using matrix_t = Eigen::Matrix<double, Dynamic, Dynamic>;
  using sparse_matrix_t = Eigen::SparseMatrix<double>;
  using sparse_solver_t = Eigen::SparseLU<sparse_matrix_t>;
  
  iterative_method(const Triangulation<N, N> &domain, double mu, double alpha,
                const BinaryMatrix<Dynamic>& incidence_matrix, const vector_t &time_locations,
                const vector_t &observations)
      : domain_(domain), mu_(mu), alpha_(alpha), incidence_matrix_(incidence_matrix),
        time_locations_(time_locations),
        dT_(time_locations[1] - time_locations[0]), observations_(observations),
        n_locs_(incidence_matrix.rows()), n_time_locs_(time_locations.size()) {
    init();
  }

  iterative_method(const Triangulation<N, N> &domain, double mu, double alpha,
                const matrix_t &locations, const vector_t &time_locations,
                const vector_t &observations)
      : domain_(domain), mu_(mu), alpha_(alpha), locations_(locations),
        time_locations_(time_locations),
        dT_(time_locations[1] - time_locations[0]), observations_(observations),
        n_locs_(locations.rows()), n_time_locs_(time_locations.size()) {
    init();
  }


  void solve(double lambda,
            LBFGS<Dynamic> opt_ = LBFGS<Dynamic>(), 
            BacktrackingLineSearch line_search = BacktrackingLineSearch(2.0,0.5,1e-4)) {
    lambda_ = lambda;
    u_ = opt_.optimize(*this, u_init_, line_search); //
    n_iter_ = opt_.n_iter();
    y_ = state(u_);
    p_ = adjoint(y_);
  }

  double J(const vector_t &y, const vector_t &u) const {

    double SSE = 0;
    for (int t = 0; t < n_time_locs_; ++t) {
      SSE += (observations_.block(n_locs_ * t, 0, n_locs_, 1) - 
                Psi_ * y.block(n_dofs_ * t, 0, n_dofs_, 1)).squaredNorm(); //get_t(y, t)
    }
    return 1. / (n_locs_ * n_time_locs_) * SSE + lambda_ * u.squaredNorm();
  }

  vector_t state(const vector_t &u) const {

    FeSpace Vh(domain_, P1<1>);
    TrialFunction uh(Vh);
    TestFunction vh(Vh);
    FeFunction y_old(Vh);

    auto reac = integral(domain_)(alpha_ * y_old * uh * vh);

    vector_t y = vector_t::Zero(n_dofs_ * n_time_locs_);
    y.block(0, 0, n_dofs_, 1) = y0_;

    for (int t = 0; t < n_time_locs_ - 1; t++) {
      y_old = y.block(n_dofs_ * t, 0, n_dofs_, 1);//get_t(y, t);
      auto R_ = reac.assemble();
      R_.makeCompressed();

      sparse_matrix_t S = 1. / dT_ * M_ + A_ + R_;
      S.prune(1e-10);
      S.makeCompressed();

      sparse_solver_t lin_solver(S);
      lin_solver.factorize(S);

      vector_t b = 1. / dT_ * M_ * y_old.coeff() + M_ *u.block(n_dofs_*(t+1), 0, n_dofs_, 1); //get_t(u, t + 1);

      y.block((t + 1) * n_dofs_, 0, n_dofs_, 1) = lin_solver.solve(b);
    }

    return y;
  }

  vector_t adjoint(const vector_t &y) const {

    FeSpace Vh(domain_, P1<1>);
    TrialFunction uh(Vh);
    TestFunction vh(Vh);
    FeFunction y_old(Vh);

    auto reac = integral(domain_)(alpha_ * y_old * uh * vh);

    vector_t p = vector_t::Zero(n_dofs_ * n_time_locs_); // p(T) == 0

    for (int t = n_time_locs_ - 1; t > 0; t--) {
      y_old = y.block(n_dofs_ * t, 0, n_dofs_, 1); 
      auto R_ = reac.assemble();
      R_.makeCompressed();

      sparse_matrix_t S = 1. / dT_ * M_ + A_ + 2. * R_;
      S.prune(1e-10);
      S.makeCompressed();

      sparse_solver_t lin_solver(S);
      lin_solver.factorize(S);

      vector_t b =
          1. / dT_ * M_ * p.block(n_dofs_ * t, 0, n_dofs_, 1) -
          1. / (n_time_locs_ * n_locs_*lambda_) * PsiTDPsi_ *y.block(n_dofs_*(t-1), 0, n_dofs_, 1) +
          1. / (n_time_locs_ * n_locs_*lambda_) * Psi_.transpose() * D_ * observations_.block(n_locs_ * (t-1), 0, n_locs_, 1);

      p.block((t - 1) * n_dofs_, 0, n_dofs_, 1) = lin_solver.solve(b);
    }

    return p;
  }

  double operator()(const vector_t &u) const {

    auto y = state(u);
    return J(y, u);
  }

  std::function<vector_t(const vector_t&)> gradient() {
    return [this](const vector_t &u) {
      vector_t y = state(u);
      vector_t p = adjoint(y);
      vector_t grad_u = Mt_ * u - Mt_ * p; 
      return grad_u;
    };
  }

  vector_t y() const { return y_; }

  vector_t u() const { return u_; }

  vector_t p() const { return p_; }

  vector_t get_t(const vector_t &v, int t) const {
    return v.block(n_dofs_ * t, 0, n_dofs_, 1);
  }

  //vector_t y(int t) const { return get_t(y(), t); }
  //vector_t u(int t) const { return get_t(u(), t); }
  //vector_t p(int t) const { return get_t(p(), t); }
  //vector_t obs(int t) const {
//     return observations_.block(n_locs_ * t, 0, n_locs_, 1);
//   }

  double J() const { 
    
    double res = (observations_ - Psi_t_ * y_).squaredNorm() + lambda_ * u_.squaredNorm();
    return res;  
  }

  int n_iter() const { return n_iter_;}

  // setters
  void set_state_initial_condition(const vector_t &y0) { y0_ = y0; }

  void set_g_initial_guess(const vector_t &u0) { u_init_ = u0; }

  // optimization algorithm custom stopping criterion
  template <typename OptimizerType> bool stop_if(OptimizerType &opt) {

    bool stop;

    vector_t y_old = state(opt.x_old);
    double loss_old = J(opt.x_old, y_old);

    vector_t y_new = state(opt.x_new);
    double loss_new = J(opt.x_new, y_new);

    stop = std::abs((loss_new - loss_old) / loss_old) < tol_;
   
    return stop;
  }

  // attribute
  const Triangulation<N, N> &domain_;

  double mu_, alpha_;
  sparse_matrix_t A_; // [A]_{ij} = \int_D mu * grad(Psi_j)*grad(Psi_i) - alpha
                       // * Psi_j * Psi_i
  sparse_matrix_t M_; // Mass

  sparse_matrix_t Psi_;
  sparse_matrix_t PsiTDPsi_;
  sparse_matrix_t D_;
  sparse_matrix_t Mt_; // Spatio-Temporal Mass matrix
  sparse_matrix_t Psi_t_;

  int n_dofs_ = 0;
  int n_locs_ = 0;
  int n_time_locs_ = 0;
  double dT_ = 0.;

  matrix_t locations_;
  vector_t time_locations_;
  vector_t observations_;
  BinaryMatrix<Dynamic> incidence_matrix_;

  double lambda_ = 1.;

  vector_t y_;
  vector_t u_;
  vector_t p_;

  vector_t y0_;
  vector_t u_init_;
  double tol_ = 1e-5;
  int n_iter_;
    
};

}

#endif // __FISHER_KPP_ITERATIVE_H__