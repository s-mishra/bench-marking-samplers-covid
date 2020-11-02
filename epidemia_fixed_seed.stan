functions {

vector reverse(vector vec) {
    int K = rows(vec);
    vector[K] rev;
    for (k in 1:K) {
        rev[k] = vec[K-k+1];
    }
    return rev;
}


  /* for multiple .stan files */
  
  /** 
   * Create group-specific block-diagonal Cholesky factor, see section 2 of
   * https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
   * @param len_theta_L An integer indicating the length of returned vector, 
   *   which lme4 denotes as m
   * @param p An integer array with the number variables on the LHS of each |
   * @param dispersion Scalar standard deviation of the errors, calles sigma by lme4
   * @param tau Vector of scale parameters whose squares are proportional to the 
   *   traces of the relative covariance matrices of the group-specific terms
   * @param scale Vector of prior scales that are multiplied by elements of tau
   * @param zeta Vector of positive parameters that are normalized into simplexes
   *   and multiplied by the trace of the covariance matrix to produce variances
   * @param rho Vector of radii in the onion method for creating Cholesky factors
   * @param z_T Vector used in the onion method for creating Cholesky factors
   * @return A vector that corresponds to theta in lme4
   */
  vector make_theta_L(int len_theta_L, int[] p, real dispersion,
                      vector tau, vector scale, vector zeta,
                      vector rho, vector z_T) {
    vector[len_theta_L] theta_L;
    int zeta_mark = 1;
    int rho_mark = 1;
    int z_T_mark = 1;
    int theta_L_mark = 1;

    // each of these is a diagonal block of the implicit Cholesky factor
    for (i in 1:size(p)) { 
      int nc = p[i];
      if (nc == 1) { // "block" is just a standard deviation
        theta_L[theta_L_mark] = tau[i] * scale[i] * dispersion;
        // unlike lme4, theta[theta_L_mark] includes the dispersion term in it
        theta_L_mark += 1;
      }
      else { // block is lower-triangular               
        matrix[nc,nc] T_i; 
        real std_dev;
        real T21;
        real trace_T_i = square(tau[i] * scale[i] * dispersion) * nc;
        vector[nc] pi = segment(zeta, zeta_mark, nc); // gamma(zeta | shape, 1)
        pi /= sum(pi);                            // thus dirichlet(pi | shape)
        
        // unlike lme4, T_i includes the dispersion term in it
        zeta_mark += nc;
        std_dev = sqrt(pi[1] * trace_T_i);
        T_i[1,1] = std_dev;
        
        // Put a correlation into T_i[2,1] and scale by std_dev
        std_dev = sqrt(pi[2] * trace_T_i);
        T21 = 2.0 * rho[rho_mark] - 1.0;
        rho_mark += 1;
        T_i[2,2] = std_dev * sqrt(1.0 - square(T21));
        T_i[2,1] = std_dev * T21;
        
        for (r in 2:(nc - 1)) { // scaled onion method to fill T_i
          int rp1 = r + 1;
          vector[r] T_row = segment(z_T, z_T_mark, r);
          real scale_factor = sqrt(rho[rho_mark] / dot_self(T_row)) * std_dev;
          z_T_mark += r;
          std_dev = sqrt(pi[rp1] * trace_T_i);
          for(c in 1:r) T_i[rp1,c] = T_row[c] * scale_factor;
          T_i[rp1,rp1] = sqrt(1.0 - rho[rho_mark]) * std_dev;
          rho_mark += 1;
        }
        
        // now vech T_i
        for (c in 1:nc) for (r in c:nc) {
          theta_L[theta_L_mark] = T_i[r,c];
          theta_L_mark += 1;
        }
      }
    }
    return theta_L;
  }
  
  /** 
  * Create group-specific coefficients, see section 2 of
  * https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
  *
  * @param z_b Vector whose elements are iid normal(0,sigma) a priori
  * @param theta Vector with covariance parameters as defined in lme4
  * @param p An integer array with the number variables on the LHS of each |
  * @param l An integer array with the number of levels for the factor(s) on 
  *   the RHS of each |
  * @return A vector of group-specific coefficients
  */
  vector make_b(vector z_b, vector theta_L, int[] p, int[] l) {
    vector[rows(z_b)] b;
    int b_mark = 1;
    int theta_L_mark = 1;
    for (i in 1:size(p)) {
      int nc = p[i];
      if (nc == 1) {
        real theta_L_start = theta_L[theta_L_mark];
        for (s in b_mark:(b_mark + l[i] - 1)) 
          b[s] = theta_L_start * z_b[s];
        b_mark += l[i];
        theta_L_mark += 1;
      }
      else {
        matrix[nc,nc] T_i = rep_matrix(0, nc, nc);
        for (c in 1:nc) {
          T_i[c,c] = theta_L[theta_L_mark];
          theta_L_mark += 1;
          for(r in (c+1):nc) {
            T_i[r,c] = theta_L[theta_L_mark];
            theta_L_mark += 1;
          }
        }
        for (j in 1:l[i]) {
          vector[nc] temp = T_i * segment(z_b, b_mark, nc);
          b_mark -= 1;
          for (s in 1:nc) b[b_mark + s] = temp[s];
          b_mark += nc + 1;
        }
      }
    }
    return b;
  }

  /** 
   * Prior on group-specific parameters
   *
   * @param z_b A vector of primitive coefficients
   * @param z_T A vector of primitives for the unit vectors in the onion method
   * @param rho A vector radii for the onion method
   * @param zeta A vector of primitives for the simplexes
   * @param tau A vector of scale parameters
   * @param regularization A real array of LKJ hyperparameters
   * @param delta A real array of concentration paramters
   * @param shape A vector of shape parameters
   * @param t An integer indicating the number of group-specific terms
   * @param p An integer array with the number variables on the LHS of each |
   * @return target()
   */
  real decov_lp(vector z_b, vector z_T, vector rho, vector zeta, vector tau,
                real[] regularization, real[] delta, vector shape,
                int t, int[] p) {
    int pos_reg = 1;
    int pos_rho = 1;
    target += normal_lpdf(z_b | 0, 1);
    target += normal_lpdf(z_T | 0, 1);
    for (i in 1:t) if (p[i] > 1) {
      vector[p[i] - 1] shape1;
      vector[p[i] - 1] shape2;
      real nu = regularization[pos_reg] + 0.5 * (p[i] - 2);
      pos_reg += 1;
      shape1[1] = nu;
      shape2[1] = nu;
      for (j in 2:(p[i]-1)) {
        nu -= 0.5;
        shape1[j] = 0.5 * j;
        shape2[j] = nu;
      }
      target += beta_lpdf(rho[pos_rho:(pos_rho + p[i] - 2)] | shape1, shape2);
      pos_rho += p[i] - 1;
    }
    target += gamma_lpdf(zeta | delta, 1);
    target += gamma_lpdf(tau  | shape, 1);
    return target();
  }
  
  /**
   * Hierarchical shrinkage parameterization
   *
   * @param z_beta A vector of primitive coefficients
   * @param global A real array of positive numbers
   * @param local A vector array of positive numbers
   * @param global_prior_scale A positive real number
   * @param error_scale 1 or sigma in the Gaussian case
   * @param c2 A positive real number
   * @return A vector of coefficientes
   */
  vector hs_prior(vector z_beta, real[] global, vector[] local, 
                  real global_prior_scale, real error_scale, real c2) {
    int K = rows(z_beta);
    vector[K] lambda = local[1] .* sqrt(local[2]);
    real tau = global[1] * sqrt(global[2]) * global_prior_scale * error_scale;
    vector[K] lambda2 = square(lambda);
    vector[K] lambda_tilde = sqrt( c2 * lambda2 ./ (c2 + square(tau) * lambda2) );
    return z_beta .* lambda_tilde * tau;
  }

  /** 
   * Hierarchical shrinkage plus parameterization
   *
   * @param z_beta A vector of primitive coefficients
   * @param global A real array of positive numbers
   * @param local A vector array of positive numbers
   * @param global_prior_scale A positive real number
   * @param error_scale 1 or sigma in the Gaussian case
   * @param c2 A positive real number
   * @return A vector of coefficientes
   */
  vector hsplus_prior(vector z_beta, real[] global, vector[] local, 
                      real global_prior_scale, real error_scale, real c2) {
    int K = rows(z_beta);
    vector[K] lambda = local[1] .* sqrt(local[2]);
    vector[K] eta = local[3] .* sqrt(local[4]);
    real tau = global[1] * sqrt(global[2]) * global_prior_scale * error_scale;
    vector[K] lambda_eta2 = square(lambda .* eta);
    vector[K] lambda_tilde = sqrt( c2 * lambda_eta2 ./ 
                                 ( c2 + square(tau) * lambda_eta2) );
    return z_beta .* lambda_tilde * tau;
  }
  
  /** 
   * Cornish-Fisher expansion for standard normal to Student t
   *
   * See result 26.7.5 of
   * http://people.math.sfu.ca/~cbm/aands/page_949.htm
   *
   * @param z A scalar distributed standard normal
   * @param df A scalar degrees of freedom
   * @return An (approximate) Student t variate with df degrees of freedom
   */
  real CFt(real z, real df) {
    real z2 = square(z);
    real z3 = z2 * z;
    real z5 = z2 * z3;
    real z7 = z2 * z5;
    real z9 = z2 * z7;
    real df2 = square(df);
    real df3 = df2 * df;
    real df4 = df2 * df2;
    return z + (z3 + z) / (4 * df) + (5 * z5 + 16 * z3 + 3 * z) / (96 * df2)
           + (3 * z7 + 19 * z5 + 17 * z3 - 15 * z) / (384 * df3)
           + (79 * z9 + 776 * z7 + 1482 * z5 - 1920 * z3 - 945 * z) / (92160 * df4);
  }

  /** 
   * Return two-dimensional array of group membership
   *
   * @param N An integer indicating the number of observations
   * @param t An integer indicating the number of grouping variables
   * @param v An integer array with the indices of group membership
   * @return An two-dimensional integer array of group membership
   */
  int[,] make_V(int N, int t, int[] v) {
    int V[t,N];
    int pos = 1;
    if (t > 0) for (j in 1:N) for (i in 1:t) {
      V[i,j] = v[pos] + 1;
      pos += 1;
    }
    return V;
  }

  /** 
  * faster version of csr_matrix_times_vector
  * declared here and defined in C++
  *
  * @param m Integer number of rows
  * @param n Integer number of columns
  * @param w Vector (see reference manual)
  * @param v Integer array (see reference manual)
  * @param u Integer array (see reference manual)
  * @param b Vector that is multiplied from the left by the CSR matrix
  * @return A vector that is the product of the CSR matrix and b
  */
/*  vector csr_matrix_times_vector2(int m, int n, vector w, 
                                int[] v, int[] u, vector b);*/

  /**
   * Calculate lower bound on intercept
   *
   * @param family Integer family code
   *   1 = gaussian
   *   2 = gamma
   *   3 = inv-gaussian
   *   4 = beta
   *   5 = binomial
   *   6 = poisson
   *   7 = neg-binom
   *   8 = poisson w/ gamma noise (not currently used but in count.stan)
   * @param link Integer link code
   * @return real lower bound
   */
  real make_lower(int family, int link) {
    if (family == 1) return negative_infinity(); // Gaussian
    if (family <= 3) { // Gamma or inverse Gaussian
      if (link == 2) return negative_infinity(); // log
      return 0;
    }
    return negative_infinity();
  }

  /**
   * Calculate upper bound on intercept
   *
   * @param family Integer family code (see make_lower above for codes)
   * @param link Integer link code
   * @return real upper bound
   */
  real make_upper(int family, int link) {
    if (family == 4 && link == 5) return 0;
    return positive_infinity();
  }
  

  /** 
   * Apply inverse link function to linear predictor
   *
   * @param eta Linear predictor vector
   * @param link An integer indicating the link function
   * @return A vector, i.e. inverse-link(eta)
   */
  vector linkinv_gauss(vector eta, int link) {
    if (link == 1)      return eta;
    else if (link == 2) return exp(eta); 
    else if (link == 3) return inv(eta);
    else reject("Invalid link");
    return eta; // never reached
  }

  /** 
  * Pointwise (pw) log-likelihood vector
  *
  * @param y A vector corresponding to the outcome variable.
  * @param link An integer indicating the link function
  * @return A vector
  */
  vector pw_gauss(vector y, vector eta, real sigma, int link) {
    return -0.5 * log(6.283185307179586232 * sigma) - 
            0.5 * square((y - linkinv_gauss(eta, link)) / sigma);
  }

  /** 
  * Apply inverse link function to linear predictor
  *
  * @param eta Linear predictor vector
  * @param link An integer indicating the link function
  * @return A vector, i.e. inverse-link(eta)
  */
  vector linkinv_gamma(vector eta, int link) {
    if (link == 1)      return eta;
    else if (link == 2) return exp(eta);
    else if (link == 3) return inv(eta);
    else reject("Invalid link");
    return eta; // never reached
  }

  /** 
  * Pointwise (pw) log-likelihood vector
  *
  * @param y A vector corresponding to the outcome variable.
  * @param eta A vector of linear predictors
  * @param shape A real number for the shape parameter
  * @param link An integer indicating the link function
  * @param sum_log_y A scalar equal to the sum of log(y)
  * @return A scalar log-likelihood
  */
  real GammaReg(vector y, vector eta, real shape, 
                int link, real sum_log_y) {
    real ret = rows(y) * (shape * log(shape) - lgamma(shape)) +
               (shape - 1) * sum_log_y;
    if (link == 2)      // link is log
      ret -= shape * sum(eta) + shape * sum(y ./ exp(eta));
    else if (link == 1) // link is identity
      ret -= shape * sum(log(eta)) + shape * sum(y ./ eta);
    else if (link == 3) // link is inverse
      ret += shape * sum(log(eta)) - shape * dot_product(eta, y);
    else reject("Invalid link");
    return ret;
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector
  *
  * @param y A vector corresponding to the outcome variable.
  * @param shape A real number for the shape parameter
  * @param link An integer indicating the link function
  * @return A vector
  */
  vector pw_gamma(vector y, vector eta, real shape, int link) {
    int N = rows(eta);
    vector[N] ll;
    if (link == 3) { // link = inverse
      for (n in 1:N) {
        ll[n] = gamma_lpdf(y[n] | shape, shape * eta[n]);
      }
    }
    else if (link == 2) { // link = log
      for (n in 1:N) {
        ll[n] = gamma_lpdf(y[n] | shape, shape / exp(eta[n]));
      }
    }
    else if (link == 1) { // link = identity
      for (n in 1:N) {
        ll[n] = gamma_lpdf(y[n] | shape, shape / eta[n]);
      }
    }
    else reject("Invalid link");
    return ll;
  }

  /** 
  * Apply inverse link function to linear predictor
  *
  * @param eta Linear predictor vector
  * @param link An integer indicating the link function
  * @return A vector, i.e. inverse-link(eta)
  */
  vector linkinv_inv_gaussian(vector eta, int link) {
    if (link == 1)      return eta;
    else if (link == 2) return exp(eta);
    else if (link == 3) return inv(eta);
    else if (link == 4) return inv_sqrt(eta);
    else reject("Invalid link");
    return eta; // never reached
  }

  /** 
  * inverse Gaussian log-PDF
  *
  * @param y The vector of outcomes
  * @param mu The vector of conditional means
  * @param lambda A positive scalar dispersion parameter
  * @param sum_log_y A scalar equal to the sum of log(y)
  * @param sqrt_y A vector equal to sqrt(y)
  * @return A scalar
  */
  real inv_gaussian(vector y, vector mu, real lambda, 
                    real sum_log_y, vector sqrt_y) {
    return 0.5 * rows(y) * log(lambda / 6.283185307179586232) - 
      1.5 * sum_log_y - 
      0.5 * lambda * dot_self( (y - mu) ./ (mu .* sqrt_y) );
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector
  *
  * @param y A vector corresponding to the outcome variable.
  * @param eta The linear predictors
  * @param lamba A positive scalar dispersion parameter
  * @param link An integer indicating the link function
  * @param log_y A precalculated vector of the log of y
  * @param sqrt_y A precalculated vector of the square root of y
  * @return A vector of log-likelihoods
  */
  vector pw_inv_gaussian(vector y, vector eta, real lambda, 
                         int link, vector log_y, vector sqrt_y) {
    vector[rows(y)] mu = linkinv_inv_gaussian(eta, link); // link checked
    return -0.5 * lambda * square( (y - mu) ./ (mu .* sqrt_y) ) +
            0.5 * log(lambda / 6.283185307179586232) - 1.5 * log_y;
  }
  
  /** 
  * PRNG for the inverse Gaussian distribution
  *
  * Algorithm from wikipedia 
  *
  * @param mu The expectation
  * @param lambda The dispersion
  * @return A draw from the inverse Gaussian distribution
  */
  real inv_gaussian_rng(real mu, real lambda) {
    real mu2 = square(mu);
    real z = uniform_rng(0,1);
    real y = square(normal_rng(0,1));
    real x = mu + ( mu2 * y - mu * sqrt(4 * mu * lambda * y + mu2 * square(y)) )
           / (2 * lambda);
    if (z <= (mu / (mu + x))) return x;
    else return mu2 / x;
  }
  
  /** 
  * Apply inverse link function to linear predictor for beta models
  *
  * @param eta Linear predictor vector
  * @param link An integer indicating the link function
  * @return A vector, i.e. inverse-link(eta)
  */
  vector linkinv_beta(vector eta, int link) {
    if (link == 1) return inv_logit(eta);  // logit
    else if (link == 2) return Phi(eta);   // probit
    else if (link == 3) return inv_cloglog(eta);  // cloglog
    else if (link == 4) return 0.5 + atan(eta) / pi(); // cauchy
    else if (link == 5) return exp(eta); // log 
    else if (link == 6) return 1 - inv_cloglog(-eta); // loglog
    else reject("invalid link");
    return eta; // never reached
  }
  
  /** 
  * Apply inverse link function to linear predictor for dispersion for beta models
  *
  * @param eta Linear predictor vector
  * @param link An integer indicating the link function
  * @return A vector, i.e. inverse-link(eta)
  */
  vector linkinv_beta_z(vector eta, int link) {
    if (link == 1) return exp(eta);         // log
    else if (link == 2) return eta;         // identity
    else if (link == 3) return square(eta); // sqrt
    else reject("Invalid link");
    return eta; // never reached
  }
  
  /** 
  * Pointwise (pw) log-likelihood vector for beta models
  *
  * @param y The vector of outcomes
  * @param eta The linear predictors
  * @param dispersion Positive dispersion parameter
  * @param link An integer indicating the link function
  * @return A vector of log-likelihoods
  */
  vector pw_beta(vector y, vector eta, real dispersion, int link) {
    vector[rows(y)] ll;
    vector[rows(y)] mu = linkinv_beta(eta, link); // link checked
    for (n in 1:rows(y)) {
      ll[n] = beta_lpdf(y[n] | mu[n] * dispersion, (1 - mu[n]) * dispersion);
    }
    return ll;
  }

  /** 
  * Pointwise (pw) log-likelihood vector for beta models with z variables
  *
  * @param y The vector of outcomes
  * @param eta The linear predictors (for y)
  * @param eta_z The linear predictors (for dispersion)
  * @param link An integer indicating the link function passed to linkinv_beta
  * @param link_phi An integer indicating the link function passed to linkinv_beta_z
  * @return A vector of log-likelihoods
  */
  vector pw_beta_z(vector y, vector eta, vector eta_z, int link, int link_phi) {
    vector[rows(y)] ll;
    vector[rows(y)] mu = linkinv_beta(eta, link); // link checked
    vector[rows(y)] mu_z = linkinv_beta_z(eta_z, link_phi); // link checked
    for (n in 1:rows(y)) {
      ll[n] = beta_lpdf(y[n] | mu[n] * mu_z[n], (1-mu[n]) * mu_z[n]);
    }
    return ll;
  }
  

 /** 
   * Apply inverse link function to linear predictor.
   * Adapted from rstanarm
   *
   * @param eta Linear predictor vector
   * @param link An integer indicating the link function
   * @return A vector, i.e. inverse-link(eta)
*/
  vector linkinv(vector eta, int link) {
    if (link == 1)      return(inv_logit(eta)); // logit
    else if (link == 2) return(Phi(eta)); // probit
    else if (link == 3) return(atan(eta) / pi() + 0.5);  // cauchit
    else if (link == 4) return(inv_cloglog(eta)); // cloglog
    else if (link == 5) return(eta); //identity
    else reject("Invalid link");
    return eta; // never reached
  }

  vector test_csr_matrix_times_vector(int m, int n, vector w,
                                      int[] v, int[] u, vector b) {
    return csr_matrix_times_vector(m, n, w, v, u, b);
  }
}

data {


int<lower=1> M; // number of countries
int<lower=1> N0; // number of time points for which to impute infections
int<lower=1> starts[M]; // the start index of each group
int<lower=1> NC[M]; // days of observed data for each group.
int<lower=1> N; // sum of NC
int<lower=1> N2; // total period for the simulation
int<lower=1> NS; // maximum number of simulation days for any given group



int<lower=0> N_obs; // total size of the observation vector
int obs_group[N_obs]; // group (1 to M) to which each observation belongs
int obs_type[N_obs]; // type of observation (1 to r). 
int obs_date[N_obs]; // observation date (1 to N2)
int y[M]; // fixed seeding for each country

int<lower=0> R; // number of different observation types
int<lower=0> oN[10];  // number of each observation type
int<lower=1> pvecs_len[R]; // maximum lag for each i2o distribution
vector<lower=0>[NS] pvecs[R]; // the 'i2o' for each type of observation

int<lower=0, upper=1> has_offset[R];
vector[N_obs] offset_;

// family for each observation type
int<lower=1,upper=3> ofamily[R]; //1:poisson 2:neg_binom 3:quasi-poisson
int<lower=1,upper=5> olink[R]; //1:log 2:probit 3:cauchit 4:cloglog 5:identity
  
// data for auxiliary parameters
int<lower=0> num_oaux; // total number aux params
int<lower=0, upper=num_oaux> has_oaux[R];
int obs[N_obs]; // vector of observations


// data relating to model matrices for each observation
int<lower=0> oK[10];  // number of predictors for each observation type
int<lower=0> K_all; // sum of the above
int<lower=0> num_ointercepts; // total intercept params
int<lower=0, upper=num_ointercepts> has_ointercept[R]; // 0 means no, otherwise gives index

// model matrices (maximum of 10 types)
// not pretty, but hopefully more efficient with algorithmic diff
vector[K_all] oxbar;
matrix[oN[1],oK[1]] oX1; 
matrix[oN[2],oK[2]] oX2;
matrix[oN[3],oK[3]] oX3;
matrix[oN[4],oK[4]] oX4;
matrix[oN[5],oK[5]] oX5;
matrix[oN[6],oK[6]] oX6;
matrix[oN[7],oK[7]] oX7;
matrix[oN[8],oK[8]] oX8;
matrix[oN[9],oK[9]] oX9;
matrix[oN[10],oK[10]] oX10;

// priors for model parameters
real<lower=0> r0; // the prior expected value for r0
vector<lower=0>[M] pop;
int<lower=1> si_len;
simplex[NS] si; // fixed serial interval using empirical data
int<lower=0, upper=1> pop_adjust; // whether to perform population adjustment or not





// dimensions
int<lower=0> K;  // number of predictors
// data
vector[K] xbar;      // predictor means
matrix[N,K] X;       // centered predictor matrix in the dense case



// flag indicating whether to draw from the prior
int<lower=0,upper=1> prior_PD;  // 1 = yes

// intercept
int<lower=0,upper=1> has_intercept;  // 1 = yes

// prior family: 0 = none, 1 = normal, 2 = student_t, 3 = hs, 4 = hs_plus, 
//   5 = laplace, 6 = lasso, 7 = product_normal
int<lower=0,upper=8> prior_dist;
int<lower=0,upper=2> prior_dist_for_intercept[has_intercept];

int<lower=0> ac_nterms; // total number of calls to rw() in the formula
int<lower=0> ac_nproc; // total number of autocorrelation processes
int<lower=0> ac_q; // total number of time periods (sum of time periods for each process)
int<lower=0> ac_nnz; // total number of non-zero entries
int<lower=0> ac_ntime[ac_nproc]; // number of time periods for each process
int<lower=-1, upper=ac_q-1> ac_v[ac_nnz]; // column indices from rstan::extract_sparse_matrix, -1 corresponds to no ac_term


// hyperparameter values are set to 0 if there is no prior
vector<lower=0>[K] prior_scale;
real<lower=0> prior_scale_for_intercept[has_intercept];
vector[K] prior_mean;
vector<lower=0>[K] prior_shape;
vector[K] prior_shift;
real prior_mean_for_intercept[has_intercept];
vector<lower=0>[K] prior_df;
real<lower=0> prior_df_for_intercept[has_intercept];
real<lower=0> global_prior_df;     // for hs priors only
real<lower=0> global_prior_scale;  // for hs priors only
real<lower=0> slab_df;     // for hs prior only
real<lower=0> slab_scale;  // for hs prior only
int<lower=2> num_normals[prior_dist == 7 ? K : 0];

// additional hyperparameters for coefficients in obs regressions
vector[K_all] prior_omean;
vector<lower=0>[K_all] prior_oscale;
vector[num_ointercepts] prior_mean_for_ointercept;
vector<lower=0>[num_ointercepts] prior_scale_for_ointercept;

// and also for auxiliary variables
int<lower=0, upper=3> prior_dist_for_oaux[num_oaux];
vector[num_oaux] prior_mean_for_oaux;
vector<lower=0>[num_oaux] prior_scale_for_oaux;
vector<lower=0>[num_oaux] prior_df_for_oaux;

real<lower=0> prior_scale_for_tau;
vector<lower=0>[ac_nproc] ac_prior_scales; // prior scale for hyperparameter for each walk.

  // glmer stuff, see table 3 of
  // https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
  int<lower=0> t;               // num. terms (maybe 0) with a | in the glmer formula
  int<lower=1> p[t];            // num. variables on the LHS of each |
  int<lower=1> l[t];            // num. levels for the factor(s) on the RHS of each |
  int<lower=0> q;               // conceptually equals \sum_{i=1}^t p_i \times l_i
  int<lower=0> len_theta_L;     // length of the theta_L vector

  // hyperparameters for glmer stuff; if t > 0 priors are mandatory
  vector<lower=0>[t] shape; 
  vector<lower=0>[t] scale;
  int<lower=0> len_concentration;
  real<lower=0> concentration[len_concentration];
  int<lower=0> len_regularization;
  real<lower=0> regularization[len_regularization];

  int<lower=0> num_non_zero;  // number of non-zero elements in the Z matrix
  vector[num_non_zero] w;     // non-zero elements in the implicit Z matrix
  int<lower=0, upper=q-1> v[num_non_zero];               // column indices for w
  int<lower=0, upper=rows(w) + 1> u[t > 0 ? N + 1 : 0];  // where the non-zeros start in each row
  int<lower=0,upper=1> special_case;                     // is the only term (1|group)
}

transformed data {
  real aux = not_a_number();
  int<lower=1> V[special_case ? t : 0, N] = make_V(N, special_case ? t : 0, v);
  int<lower=0> ac_V[ac_nterms, N] = make_V(N, ac_nterms, ac_v);


simplex[NS] si_rev = reverse(si); 
vector<lower=0>[NS] pvecs_rev[R]; 


  int<lower=0> len_z_T = 0;
  int<lower=0> len_var_group = sum(p) * (t > 0);
  int<lower=0> len_rho = sum(p) - t;
  int<lower=1> pos = 1;
  real<lower=0> delta[len_concentration];
  int<lower=0> hs;
  if (prior_dist <= 2) hs = 0;
  else if (prior_dist == 3) hs = 2;
  else if (prior_dist == 4) hs = 4;
  else hs = 0;
  
  for (i in 1:t) {
    if (p[i] > 1) {
      for (j in 1:p[i]) {
        delta[pos] = concentration[j];
        pos += 1;
      }
    }
    for (j in 3:p[i]) len_z_T += p[i] - 1;
  }

for(r in 1:R)
      pvecs_rev[r] = reverse(pvecs[r]);

}

parameters {
  vector[num_ointercepts] ogamma;
  real gamma[has_intercept];
  vector<lower=0>[num_oaux] oaux_raw;

vector<upper=(prior_dist == 8 ? 0 : positive_infinity())>[prior_dist == 7 ? sum(num_normals) : K] z_beta;
real<lower=0> global[hs];
vector<lower=0>[K] local[hs];
real<lower=0> caux[hs > 0];
vector<lower=0>[K] mix[prior_dist == 5 || prior_dist == 6];
real<lower=0> one_over_lambda[prior_dist == 6];
vector[q] z_b;
vector[len_z_T] z_T;
vector<lower=0,upper=1>[len_rho] rho;
vector<lower=0>[len_concentration] zeta;
vector<lower=0>[t] tau;


vector<lower=0>[ac_nproc] ac_scale_raw; // standard normal version of scale hyperparameter
vector[ac_q] ac_noise; // noise terms for each walk at each time period

vector[K_all] oz_beta;
}

transformed parameters {
  vector[N_obs] oeta;
  vector[N_obs] E_obs; // expected values of the observations 
  vector[N] eta;  // linear predictor
  vector<lower=0>[num_oaux] oaux = oaux_raw;


matrix<lower=0>[N2, M] Rt_unadj = rep_matrix(0,N2,M);
matrix<lower=0>[N2, M] Rt = rep_matrix(0,N2,M);
matrix<lower=0>[N2, M] infections = rep_matrix(0,N2,M);
matrix<lower=0>[N2, M] infectiousness = rep_matrix(0,N2,M);

vector[ac_nproc] ac_scale = ac_scale_raw .* ac_prior_scales;
vector[ac_q] ac_beta;


vector[K_all] obeta = oz_beta .* prior_oscale + prior_omean;

  vector[K] beta;
  vector[q] b;
  vector[len_theta_L] theta_L;
  if      (prior_dist == 0) beta = z_beta;
  else if (prior_dist == 1) beta = z_beta .* prior_scale + prior_mean;
  else if (prior_dist == 2) for (k in 1:K) {
    beta[k] = CFt(z_beta[k], prior_df[k]) * prior_scale[k] + prior_mean[k];
  }
  else if (prior_dist == 3) {
    real c2 = square(slab_scale) * caux[1];
      beta = hs_prior(z_beta, global, local, global_prior_scale, 1, c2);
  }
  else if (prior_dist == 4) {
    real c2 = square(slab_scale) * caux[1];
      beta = hsplus_prior(z_beta, global, local, global_prior_scale, 1, c2);
  }
  else if (prior_dist == 5) // laplace
    beta = prior_mean + prior_scale .* sqrt(2 * mix[1]) .* z_beta;
  else if (prior_dist == 6) // lasso
    beta = prior_mean + one_over_lambda[1] * prior_scale .* sqrt(2 * mix[1]) .* z_beta;
  else if (prior_dist == 7) { // product_normal
    int z_pos = 1;
    for (k in 1:K) {
      beta[k] = z_beta[z_pos];
      z_pos += 1;
      for (n in 2:num_normals[k]) {
        beta[k] *= z_beta[z_pos];
        z_pos += 1;
      }
      beta[k] *= prior_scale[k] ^ num_normals[k];
      beta[k] += prior_mean[k];
    }
  }
  else if(prior_dist == 8) beta = z_beta + prior_shift; // shifted gamma

// transform auxiliary parameters
  for (i in 1:num_oaux) {
    if (prior_dist_for_oaux[i] > 0) {
      if (prior_scale_for_oaux[i] > 0) {
        oaux[i] *= prior_scale_for_oaux[i];
      }
      if (prior_dist_for_oaux[i] <= 2) {
        oaux[i] += prior_mean_for_oaux[i];
      }
    }
  }

  {
    int i = 1;
    for (proc in 1:ac_nproc) { // this treats ac terms as random walks for now (to be extended to AR(p))
        ac_beta[i:(i+ac_ntime[proc]-1)] = cumulative_sum(ac_noise[i:(i+ac_ntime[proc]-1)]);
        i += ac_ntime[proc];
    }
  }
  

  if (K > 0) {
     eta = X * beta;
  }
  else eta = rep_vector(0.0, N);

  if (t > 0) {
    if (special_case == 1) {
      int start = 1;
      theta_L = scale .* tau;
      if (t == 1) b = theta_L[1] * z_b;
      else for (i in 1:t) {
        int end = start + l[i] - 1;
        b[start:end] = theta_L[i] * z_b[start:end];
        start = end + 1;
      }
    }
    else {
      theta_L = make_theta_L(len_theta_L, p,
                             1.0, tau, scale, zeta, rho, z_T);
      b = make_b(z_b, theta_L, p, l);
    }
  }

  if (t > 0) {

if (special_case) for (i in 1:t) eta += b[V[i]];
else eta += csr_matrix_times_vector(N, q, w, v, u, b);
  }
  if (has_intercept == 1) {
    eta += gamma[1];
  }
  else {

  // correction to eta if model has no intercept (because X is centered)
  eta += dot_product(xbar, beta); 
  }

  if (ac_nterms > 0) {

for (i in 1:ac_nterms) {
    for (j in 1:N) {
        if (ac_V[i,j] > 0) { # if 0 then doesn't have RW component
            eta[j] += ac_beta[ac_V[i,j]];
        }
    }
}
  }
  


{
    int npos = 1;
    int kpos = 1;
    int i = 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            oeta[npos:(npos+oN[i]-1)] = oX1 * segment(obeta, kpos, oK[i]);
        else {
            oeta[npos:(npos+oN[i]-1)] = rep_vector(0.0, oN[i]);
        }
        // add intercept
        if (has_ointercept[i] > 0)
            oeta[npos:(npos+oN[i]-1)] += ogamma[has_ointercept[i]];
        else if (oK[i] > 0) {
            oeta[npos:(npos+oN[i]-1)] += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]));
        }
    }

    npos += oN[i];
    kpos += oK[i];
    i += 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            oeta[npos:(npos+oN[i]-1)] = oX2 * segment(obeta, kpos, oK[i]);
        else {
            oeta[npos:(npos+oN[i]-1)] = rep_vector(0.0, oN[i]);
        }
        // add intercept
        if (has_ointercept[i] > 0)
            oeta[npos:(npos+oN[i]-1)] += ogamma[has_ointercept[i]];
        else if (oK[i] > 0) {
            oeta[npos:(npos+oN[i]-1)] += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]));
        }
    }

    npos += oN[i];
    kpos += oK[i];
    i += 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            oeta[npos:(npos+oN[i]-1)] = oX3 * segment(obeta, kpos, oK[i]);
        else {
            oeta[npos:(npos+oN[i]-1)] = rep_vector(0.0, oN[i]);
        }
        // add intercept
        if (has_ointercept[i] > 0)
            oeta[npos:(npos+oN[i]-1)] += ogamma[has_ointercept[i]];
        else if (oK[i] > 0) {
            oeta[npos:(npos+oN[i]-1)] += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]));
        }
    }

    npos += oN[i];
    kpos += oK[i];
    i += 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            oeta[npos:(npos+oN[i]-1)] = oX4 * segment(obeta, kpos, oK[i]);
        else {
            oeta[npos:(npos+oN[i]-1)] = rep_vector(0.0, oN[i]);
        }
        // add intercept
        if (has_ointercept[i] > 0)
            oeta[npos:(npos+oN[i]-1)] += ogamma[has_ointercept[i]];
        else if (oK[i] > 0) {
            oeta[npos:(npos+oN[i]-1)] += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]));
        }
    }

    npos += oN[i];
    kpos += oK[i];
    i += 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            oeta[npos:(npos+oN[i]-1)] = oX5 * segment(obeta, kpos, oK[i]);
        else {
            oeta[npos:(npos+oN[i]-1)] = rep_vector(0.0, oN[i]);
        }
        // add intercept
        if (has_ointercept[i] > 0)
            oeta[npos:(npos+oN[i]-1)] += ogamma[has_ointercept[i]];
        else if (oK[i] > 0) {
            oeta[npos:(npos+oN[i]-1)] += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]));
        }
    }

    npos += oN[i];
    kpos += oK[i];
    i += 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            oeta[npos:(npos+oN[i]-1)] = oX6 * segment(obeta, kpos, oK[i]);
        else {
            oeta[npos:(npos+oN[i]-1)] = rep_vector(0.0, oN[i]);
        }
        // add intercept
        if (has_ointercept[i] > 0)
            oeta[npos:(npos+oN[i]-1)] += ogamma[has_ointercept[i]];
        else if (oK[i] > 0) {
            oeta[npos:(npos+oN[i]-1)] += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]));
        }
    }

    npos += oN[i];
    kpos += oK[i];
    i += 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            oeta[npos:(npos+oN[i]-1)] = oX7 * segment(obeta, kpos, oK[i]);
        else {
            oeta[npos:(npos+oN[i]-1)] = rep_vector(0.0, oN[i]);
        }
        // add intercept
        if (has_ointercept[i] > 0)
            oeta[npos:(npos+oN[i]-1)] += ogamma[has_ointercept[i]];
        else if (oK[i] > 0) {
            oeta[npos:(npos+oN[i]-1)] += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]));
        }
    }

    npos += oN[i];
    kpos += oK[i];
    i += 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            oeta[npos:(npos+oN[i]-1)] = oX8 * segment(obeta, kpos, oK[i]);
        else {
            oeta[npos:(npos+oN[i]-1)] = rep_vector(0.0, oN[i]);
        }
        // add intercept
        if (has_ointercept[i] > 0)
            oeta[npos:(npos+oN[i]-1)] += ogamma[has_ointercept[i]];
        else if (oK[i] > 0) {
            oeta[npos:(npos+oN[i]-1)] += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]));
        }
    }

    npos += oN[i];
    kpos += oK[i];
    i += 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            oeta[npos:(npos+oN[i]-1)] = oX9 * segment(obeta, kpos, oK[i]);
        else {
            oeta[npos:(npos+oN[i]-1)] = rep_vector(0.0, oN[i]);
        }
        // add intercept
        if (has_ointercept[i] > 0)
            oeta[npos:(npos+oN[i]-1)] += ogamma[has_ointercept[i]];
        else if (oK[i] > 0) {
            oeta[npos:(npos+oN[i]-1)] += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]));
        }
    }

    npos += oN[i];
    kpos += oK[i];
    i += 1;
    if (R >= i) {
        // predictor without intercept
        if (oK[i] > 0)
            oeta[npos:(npos+oN[i]-1)] = oX10 * segment(obeta, kpos, oK[i]);
        else {
            oeta[npos:(npos+oN[i]-1)] = rep_vector(0.0, oN[i]);
        }
        // add intercept
        if (has_ointercept[i] > 0)
            oeta[npos:(npos+oN[i]-1)] += ogamma[has_ointercept[i]];
        else if (oK[i] > 0) {
            oeta[npos:(npos+oN[i]-1)] += 
                dot_product(segment(oxbar, kpos, oK[i]), 
                            segment(obeta, kpos, oK[i]));
        }
    }
}

{
    int i = 1;
    for (r in 1:R) {
        if (has_offset[r] == 1) {
            oeta[i:(i+oN[r]-1)] += segment(offset_, i, oN[r]);
        }
        i += oN[r];
    }
}


{ // predict cases over time
int idx=1;
matrix[N2,M] cumm_sum = rep_matrix(0,N2,M);
for (m in 1:M){
    // time indices for group m: start date, final seed date, final date
    int n0 = starts[m];
    int n1 = n0 + N0 - 1;
    int n2 = n0 + NC[m] - 1;
    int len;

    // impute unadjusted Rt from the linear predictor
    Rt_unadj[n0:n2,m] = exp(eta[idx:(idx+NC[m]-1)]);

    Rt[n0:n1,m] = Rt_unadj[n0:n1,m]; 
    idx += NC[m];

    infections[n0:n1,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
    cumm_sum[n0:n1,m] = cumulative_sum(infections[n0:n1,m]);

    for(i in (n1+1):n2) {
        int start = max(n0, i - si_len);
        real convolution = dot_product(sub_col(infections, start, m, i - start), tail(si_rev, i - start));
        if (pop_adjust) {
            infections[i,m] = (pop[m] - cumm_sum[i-1,m]) * (1 - exp(-Rt_unadj[i,m] * convolution / pop[m]));
            Rt[i,m] =  (pop[m] - cumm_sum[i-1,m]) * Rt_unadj[i,m] / pop[m];
        } else {
            infections[i,m] = Rt_unadj[i,m] * convolution;
            Rt[i,m] = Rt_unadj[i,m];
        }
        infectiousness[i,m] = convolution / max(si);
        cumm_sum[i,m] = cumm_sum[i-1,m] + infections[i,m];
    }
}
}


{ // apply link function
    int i = 1;
    for (r in 1:R) {
        E_obs[i:(i+oN[r]-1)] = linkinv(segment(oeta, i, oN[r]) + 1e-15, olink[r]);
        i += oN[r];
    }
}
 
{  // compute expected values of the observations
for (i in 1:N_obs) {
    int m = obs_group[i];
    int dt = obs_date[i];
    int tp = obs_type[i];
    int n0 = starts[m];
    if (dt == 1)
        E_obs[i] *= 1e-15 * infections[1,m];
    else {
        int start = max(n0, dt - pvecs_len[tp]);
        E_obs[i] *= dot_product(sub_col(infections, start, m, dt-start), tail(pvecs_rev[tp], dt-start));
    }
}
}
}

model {
  // Log-priors for coefficients
       if (prior_dist == 1) target += normal_lpdf(z_beta | 0, 1);
  else if (prior_dist == 2) target += normal_lpdf(z_beta | 0, 1); // Student t via Cornish-Fisher expansion
  else if (prior_dist == 3) { // hs
    real log_half = -0.693147180559945286;
    target += normal_lpdf(z_beta | 0, 1);
    target += normal_lpdf(local[1] | 0, 1) - log_half;
    target += inv_gamma_lpdf(local[2] | 0.5 * prior_df, 0.5 * prior_df);
    target += normal_lpdf(global[1] | 0, 1) - log_half;
    target += inv_gamma_lpdf(global[2] | 0.5 * global_prior_df, 0.5 * global_prior_df);
    target += inv_gamma_lpdf(caux | 0.5 * slab_df, 0.5 * slab_df);
  }
  else if (prior_dist == 4) { // hs+
    real log_half = -0.693147180559945286;
    target += normal_lpdf(z_beta | 0, 1);
    target += normal_lpdf(local[1] | 0, 1) - log_half;
    target += inv_gamma_lpdf(local[2] | 0.5 * prior_df, 0.5 * prior_df);
    target += normal_lpdf(local[3] | 0, 1) - log_half;
    // unorthodox useage of prior_scale as another df hyperparameter
    target += inv_gamma_lpdf(local[4] | 0.5 * prior_scale, 0.5 * prior_scale);
    target += normal_lpdf(global[1] | 0, 1) - log_half;
    target += inv_gamma_lpdf(global[2] | 0.5 * global_prior_df, 0.5 * global_prior_df);
    target += inv_gamma_lpdf(caux | 0.5 * slab_df, 0.5 * slab_df);
  }
  else if (prior_dist == 5) { // laplace
    target += normal_lpdf(z_beta | 0, 1);
    target += exponential_lpdf(mix[1] | 1);
  }
  else if (prior_dist == 6) { // lasso
    target += normal_lpdf(z_beta | 0, 1);
    target += exponential_lpdf(mix[1] | 1);
    target += chi_square_lpdf(one_over_lambda[1] | prior_df[1]);
  }
  else if (prior_dist == 7) { // product_normal
    target += normal_lpdf(z_beta | 0, 1);
  }
  else if (prior_dist == 8) { // shifted gamma
    target += gamma_lpdf(-z_beta | prior_shape, 1.0 ./ prior_scale);
  }
  /* else prior_dist is 0 and nothing is added */
  
  // Log-prior for intercept  
  if (has_intercept == 1) {
    if (prior_dist_for_intercept[1] == 1)  // normal
      target += normal_lpdf(gamma | prior_mean_for_intercept[1], prior_scale_for_intercept[1]);
    else if (prior_dist_for_intercept[1] == 2)  // student_t
      target += student_t_lpdf(gamma | prior_df_for_intercept[1], prior_mean_for_intercept[1], 
                               prior_scale_for_intercept[1]);
    /* else prior_dist is 0 and nothing is added */
  }


target += normal_lpdf(ac_scale_raw | 0, 1);

{
int i = 1;
for (proc in 1:ac_nproc) { 
    target += normal_lpdf(ac_noise[i:(i+ac_ntime[proc]-1)] | 0, ac_scale[proc]);
    i += ac_ntime[proc];
}
}

target += normal_lpdf(oz_beta | 0, 1);

// prior for the intercepts
if (num_ointercepts > 0) 
    target += normal_lpdf(ogamma | prior_mean_for_ointercept, prior_scale_for_ointercept);

// priors for auxiliary variables
for (i in 1:num_oaux) {
if (prior_dist_for_oaux[i] == 1) 
    target += normal_lpdf(oaux_raw[i] | 0, 1);
else if (prior_dist_for_oaux[i] == 2)
    target += student_t_lpdf(oaux_raw[i] | prior_df_for_oaux[i], 0, 1);
else if (prior_dist_for_oaux[i] == 3)
    target += exponential_lpdf(oaux_raw[i] | 1);
}
  if (t > 0) {
    real dummy = decov_lp(z_b, z_T, rho, zeta, tau,
                          regularization, delta, shape, t, p);
  }

  if (prior_PD == 0) {
    int i = 1;
    for (r in 1:R) {
      if (ofamily[r] == 1) { // poisson
        target += poisson_lpmf(segment(obs, i, oN[r]) | segment(E_obs, i, oN[r]) + 1e-15);
      }
      else if (ofamily[r] == 2) { // neg binom
        target += neg_binomial_2_lpmf(segment(obs, i, oN[r]) | 
        segment(E_obs, i, oN[r]) + 1e-15, oaux[has_oaux[r]]);
      } else { // quasi-poisson
        target += neg_binomial_2_lpmf(segment(obs, i, oN[r]) | 
        segment(E_obs, i, oN[r]) + 1e-15, (segment(E_obs, i, oN[r]) + 1e-15) / oaux[has_oaux[r]]);
      }
      i += oN[r];
    }
  }

}

generated quantities {
  real alpha[has_intercept];
  
  if (has_intercept == 1) {
    alpha[1] = gamma[1] - dot_product(xbar, beta);
  }
} 
