# Code to implement our proposed blocked Gibbs (BGS) algorithm to sample from the 
# posterior of HDP under a univariate conjugate Gaussian mixture model.

# source functions to sample the unnormalized global weights using rejection sampler.
source("rejectionsampler.R")

# Function to draw samples from a dirichlet distribution.
# alpha : concentration parameters
rdirich = function(alpha) {
  
  x = rgamma(n = length(alpha), shape = alpha, rate = 1)
  return(x/sum(x))
}


# Function to draw a sample from f_k(.) using Rejection Sampler
# J : population size
# log.p_k = log (\prod_j \pi_jk)
# log_u = \sum_j (log u_j)
# b0 : rate parameter of prior gamma distribution of \alpha_0
# gam : concentration parameter of base DP 
# L : max number of clusters
samp.f_k = function(N, J, log.p_k, log_u, gam, L, b0){
  
  # parameters of target density
  A = gam/L
  B = b0 - log.p_k - log_u
  
  # The target density tends to a degenerate distribution at 0 as B tends to \infty
  # Setting a threshold on B to avoid numerical issues
  if(B >= 1e+09){
    S = 1e-10
    return(S)
  }
  else {
    # parameters of the cover density
    distr = param.mixture(J = J, A = A, B = B, N = N)
    
    cond = 5
    while(cond > 0){
      
      # draw a sample from the mixture density
      S = samp.mixture(weights = distr$weights, lambda = distr$lambda, q = distr$q)
      
      # check for acceptance / rejection
      Y1 = log.f_k(t = S, J = J, A = A, B = B)
      Y2 = log.mixture_dens(t = S, m = distr$m, a = distr$a, lambda = distr$lambda, q = distr$q)
      Y = Y1 - Y2
      
      U = log(runif(1))
      
      # accounting for cases when both log f_k(S) and log cover_density(S) is -\infty, 
      # where S is the obtained sample 
      if(is.infinite(Y1) | is.na(Y1) | is.infinite(Y2) | is.na(Y2) ){
        cond = 5
      }
      else cond = U - Y
      # Accept the sample if U <=Y, else repeat the sampling step.
    }
    return(S)
  }
}


# Function to draw samples from the full conditional of phi
# Z : list containing the class assignments
# x : list containing the data
# phi.param : parameters of prior distribution of phi
update_phi = function(Z, L, x, phi.param){
  
  phi = rep(NA, L)
  
  # unlisting all observations and class labels
  x.all = unlist(x)
  Z.all = unlist(Z)
  
  # n = (n1,...,nL) where n_k = \sum_(j,i) I(Z_ji = k)
  n.group = sapply(seq_len(L), function(k) sum(Z.all == k))
  
  # xbar = (xbar_1,..., xbar_L), xbar_k =(1/n_k) * \sum_{(j,i): Z_ji = k} x_ji
  xbar = sapply(seq_len(L), function(k) mean(x.all[Z.all == k]))
  
  # xi : prior mean of phi, lambda : prior precision of phi
  # Prior distribution : phi_k follows Normal(xi, lambda^{-1})
  # Likelihood : x_ji | Z_ji=k follows Normal(phi_k, tau^{-1}), 
  # tau: precision of x_ji
  xi = phi.param[1]; lambda = phi.param[2]; tau = phi.param[3]
  
  # posterior mean and standard deviation of phi
  mean_phi = ((n.group*xbar*tau) + xi*lambda)/((n.group*tau) + lambda)
  sd_phi = 1/ sqrt((n.group*tau) + lambda)
  
  # sample phi from its prior for the unoccupied clusters
  mean_phi[n.group == 0] = xi
  sd_phi[n.group == 0] = 1/sqrt(lambda)
  
  phi = sapply(seq_len(L), function(k) rnorm(1, mean = mean_phi[k], sd = sd_phi[k]))
  
  return(phi)
}


# Function to draw samples from the full conditional of Z
# Pi : JxL matrix where Pi[j, ] denotes the weight parameter for the jth population. Here, Z_ji follows Cat(1:L, Pi[j, ])
# x : list containing the data
# phi = (phi_1, ..., phi_L) : set of parameters corresponding to the global clusters.
# phi.param : parameters of prior distribution of phi
update_Z = function(Pi, x, phi, phi.param){
  
  J = length(x); n = lengths(x)
  Z = vector(mode = "list", length = J)
  
  for(j in 1:J){
    
    L = length(phi)
    tau = phi.param[3]
    F_j = sapply(seq_len(L), function(k) dnorm(x[[j]], mean = phi[k], sd = sqrt(1/tau), log = TRUE))
    
    # Pi_j : matrix whose rows contain pi_j^{*}, 
    # where \pi_{ji,k}^{*} \propto \pi_{jk} f(x_ji | \phi_k) 
    Pi_j = matrix(log(Pi[j, ]), ncol = L, nrow = n[j], byrow = TRUE) + F_j
    Pi_j = t(apply(Pi_j, 1, function(x) exp(x - max(x))))
    Pi_j = Pi_j/rowSums(Pi_j)
    
    Z[[j]] = apply(Pi_j, 1, function(x) sample(1:L, size = 1, prob = x))
  }
  
  return(Z)
}


# Function to draw samples from the full conditional of Pi
# Z : list containing the class assignments
# t = (t_1, ..., t_L) : unnormalized global weights
update_Pi = function(Z, t){
  
  J = length(Z); L = length(t)
  Pi =  matrix(NA, nrow = J, ncol = L)
  
  for(j in 1:J){
    
    # n_j = (n_j1, ..., n_jL), n_jk = \sum_i I(Z_ji = k)
    n_j = sapply(seq_len(L), function(i) sum(Z[[j]]==i))
    
    # Pi |... follows Dirichlet (n_j + t)
    Pi.draw = rdirich(n_j + t)
    
    # setting an lower bound of 10^(-10) for \pi_jk's to avoid numerical issues
    Pi.ind = which(Pi.draw < 1e-10)
    
    if(length(Pi.ind) > 0){
      
      Pi.draw[Pi.ind] = 1e-10
      
      excess.P = sum(Pi.draw) - 1
      
      ind.max = which.max(Pi.draw)
      Pi.draw[ind.max] = Pi.draw[ind.max] - excess.P
    }
    
    Pi[j, ] = Pi.draw
  }
  return(Pi)
}


# Function to draw samples from the full conditional of t and hence Beta
# b0 : rate parameter of prior gamma distribution of \alpha_0
# gam : concentration parameter of base DP
# Pi : JxL matrix where Pi[j, ] denotes the weight parameter for the jth population.
# u = (u_1, u_2, \ldots, u_J) : vector of augmented Gamma random variables
update_Beta = function(N, Pi, gam, b0, u){
  
  log_u = sum(log(u))
  L = ncol(Pi); J = nrow(Pi); log.p = colSums(log(Pi))
  
  # log.p tends to -\infty for unoccupied clusters
  # setting a threshold of -10^10 for log.p to avoid numerical issues
  log.p = sapply(1:L, function(x) ifelse(var(Pi[,x]) == 0, -1e+10, log.p))
  
  # draw t_k from f_k(.), k = 1, 2,..., L.
  t = sapply(log.p, function(p) samp.f_k(N = N, J = J, log.p_k = p, 
                                         log_u = log_u, gam = gam, L = L, b0 = b0) )
  sum_t = sum(t)
  
  # set \beta_k = t_k / \sum_k (t_k)
  Beta = t/sum_t
  
  return(list("Beta" = Beta, "sum_t" = sum_t, "t" = t))
}

# Function to draw samples from the full conditional of augmented gamma random variables, u
# sum_t = \sum_k T_k, T_k \sim f_k
# J = population size
update_u = function(sum_t, J){
  
  # u_1, \ldots, u_J \sim Gamma(sum_t, 1) independently
  res = rgamma(n = J, shape = sum_t, rate = 1)
  
  return(res)
}


# Our proposed Blocked Gibbs Sampler 
# x : list of length J, x[[j]] contains data in the jth population
# L.max : maximum number of clusters
# gam : concentration parameter of global DP
# phi.param : parameters specifying the prior distribution of \phi
# b0 : shared concentration parameter has a gamma prior with shape = gam and rate = b0
# N : Integer to get the number of knot points (= 2N+2) for the rejection sampler, default value is 1.
# Burn.in : Burn in period of the MCMC chain
# M : number of MCMC samples required
# est.density : TRUE/FALSE indicating whether the density needs to be estimated
# y.grid : grid points for estimating the density of the populations, if est.density  = TRUE
blocked_gibbs = function(x, L.max, gam, phi.param, b0, N = 1,
                         Burn.in, M, est.density = FALSE, y.grid = NULL){
  
  J = length(x)
  
  # set initial values for running the Gibbs sampler
  Pi = matrix(1/L.max, nrow = J, ncol = L.max)
  
  xi = phi.param[1]; lambda = phi.param[2]; tau = phi.param[3]
  phi = rnorm(L.max, mean = xi, sd = 1/sqrt(lambda))
  
  t = rgamma(n = L.max, shape = gam/L.max, rate = b0)
  u = rgamma(n = J, shape = 1, rate = 1)
  
  # list to store the posterior samples
  Iterates = vector(mode = "list", length = M)
  
  for(m in 1:(M + Burn.in)){
    
    # time at the beginning
    T1 = Sys.time()
    
    # update Z
    Z = update_Z(Pi = Pi, x = x, phi = phi, phi.param = phi.param)
    
    # update phi
    phi = update_phi(Z = Z, L = L.max, x = x, phi.param = phi.param)
    
    # update Pi
    Pi = update_Pi(Z = Z, t = t)
    
    # update t and Beta 
    res = update_Beta(N = N, Pi = Pi, gam = gam, b0 = b0, u = u)
    Beta = res$Beta
    sum_t = res$sum_t
    t = res$t
    
    # update u
    u = update_u(sum_t = sum_t, J = J)
    
    # time at the end of all updates
    T2 = Sys.time()
    Tdiff =  difftime(T2, T1, units = "secs")
    
    # order the phi's in decreasing order of cluster occupancy
    n.group = sapply(seq_len(L.max), function(j) sum(unlist(Z)==j))
    dat = data.frame(phi, n.group, ind = 1:L.max)
    dat = dat[order(-dat$n.group),]
    
    f = NULL
    
    if(est.density == TRUE){
      n.grid = length(y.grid)
      
      # J x n.grid matrix to store the densities along the rows.
      f = matrix(NA, nrow = J, ncol =  n.grid)
      for(j in 1:J){
        
        # evaluate the density for each population
        f[j, ] = sapply(1:n.grid, function(ii) 
          sum( Pi[j, ] * dnorm(y.grid[ii], mean = phi, sd = sqrt(1/tau)) ))
        
      }
    }
    
    # print every 200th iteration
    if(m %% 200 == 0){
      print(paste("iteration :", m))
    }
    
    # store samples after Burn in
    if(m > Burn.in){
      Iterates[[m-Burn.in]] = list("Z" = Z, "Pi" = Pi, "phi" = phi, "Beta" = Beta, "t" = t, 
                                   "phi.ord" = dat$phi,"n.group" = dat$n.group, "indices" = dat$ind, 
                                   "u" = u, "density" = f, "time" = Tdiff)
    }
  }
  return(Iterates)
  
}
