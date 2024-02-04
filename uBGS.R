# Code to implement a truncated version of the unmarginalized blocked Gibbs (uBGS) algorithm 
# described in Section 3 of Amini et al. (2019) to sample from the posterior of HDP 
# under a univariate conjugate Gaussian mixture model.

# source functions from BGS.R since the distinct atoms have the same posterior update under both schemes
source("BGS.R")

# Function to get unique dish indices Z_ji = k_{jt_ji} for all i,j.
# k : list containing dish indices for each table of each restaurant
# t: list containing table indices for each customer in each restaurant
get_Z = function(k, t){
  
  J = length(t); n = lengths(t)
  Z = vector(mode = "list", length = J)
  
  for(j in 1:J){
    Z[[j]] = sapply(1:n[j], function(i) k[j, t[[j]][i]])
  }
  return(Z)
}


# Function to draw samples from the full conditional of local cluster labels i.e. table indices
#t.list :list containing table indices for each customer in each restaurant
# Pi : JxT0 matrix where Pi[j, ] denotes the weight parameter for the jth population. Here, t_ji follows Cat(1:T0, Pi[j, ]) apriori.
# k.mat : JxT0 matrix containing dish indices for each table in each restaurant
# x : list containing the data
# phi = (phi_1, ..., phi_L) : set of atoms corresponding to the distinct global clusters.
# phi.param : parameters of prior distribution of phi
update_t = function(k.mat, Pi, x, phi, phi.param){
  
  J = length(x); n = lengths(x); T0 = ncol(k.mat)
  t.list = vector(mode = "list", length = J)
  tau = phi.param[3]
  
  for(j in 1:J){
    
    F_j = sapply(seq_len(T0), function(t) 
      dnorm(x[[j]], mean = phi[k.mat[j, t]], sd = sqrt(1/tau), log = TRUE))
    
    Pi_j = matrix(log(Pi[j, ]), ncol = T0, nrow = n[j], byrow = TRUE) + F_j
    Pi_j = t(apply(Pi_j, 1, function(x) exp(x - max(x))))
    Pi_j = Pi_j/rowSums(Pi_j)
    
    t.list[[j]] = apply(Pi_j, 1, function(x) sample(1:T0, size = 1, prob = x))
    
  }
  
  return(t.list)
}


# Function to draw samples from the full conditional of global cluster labels i.e. dish indices
# t.list : list containing table indices for each customer in each restaurant
# T0 : truncation level for the local weights
# Beta = (Beta_1, ..., Beta_L) : global weights. Here, k_jt follows Cat(1:L, Beta) apriori.
# x : list containing the data
# phi = (phi_1, ..., phi_L) : set of atoms corresponding to the distinct global clusters.
# phi.param : parameters of prior distribution of phi
update_k = function(t.list, T0, Beta, x, phi, phi.param){
  
  J = length(x); L = length(phi)
  tau = phi.param[3]
  k.mat = matrix(NA, nrow = J, ncol = T0)
  
  for(j in 1:J){
    for(t in 1:T0){
      
      prob_k = rep(0, L)
      if(sum(t.list[[j]] == t) == 0){
        prob_k = Beta
      }
      else{
        x.tab.t = x[[j]][t.list[[j]] == t]
        
        prob_k = sapply(seq_len(L), function(k)
          log(Beta[k]) + sum(dnorm(x = x.tab.t, mean = phi[k], sd = sqrt(1/tau), log = TRUE)) )
        
        prob_k = prob_k - max(prob_k)
        prob_k = exp(prob_k)
        prob_k = prob_k/ sum(prob_k)
      }
      k.mat[j, t] = sample(1:L, size = 1, prob = prob_k)
      
    }
  }
  
  return(k.mat)
}


# Function to draw samples from the full conditional of Pi
# t.list : list containing table indices for each customer in each restaurant
# alpha0 : shared concentration parameter of the local weights
# T0 : truncation level for the local weights
update_Pi_2 = function(t.list, alpha0, T0){
  
  J = length(t.list)
  Pi =  matrix(NA, nrow = J, ncol = T0)
  Pi.prime.mat = matrix(NA, nrow = J, ncol = (T0-1))
  
  for(j in 1:J){
    
    # n_j = (n_j1, ..., n_jT0), n_jt = \sum_i I(t_ji = t)
    n_j.group = sapply(seq_len(T0), function(t) sum(t.list[[j]]==t))
    n_j = sum(n_j.group)
    
    n_j.cum = cumsum(n_j.group)
    
    # update Pi[j, ] using stick-breaking formulation
    Pi.prime = sapply(seq_len(T0-1), function(h)
      rbeta(1, shape1 = (n_j.group[h]+1), shape2 = (alpha0 + n_j - n_j.cum[h])))
    
    Pi.prime.prod = c(1, cumprod(1 - Pi.prime))
    Pi.prime[T0] = 1
    
    Pi.draw = sapply(seq_len(T0), function(i) Pi.prime[i]*Pi.prime.prod[i])
    
    
    # setting an lower bound of 10^(-10) for \pi_jk's to avoid numerical issues
    Pi.ind = which(Pi.draw < 1e-10)
    
    if(length(Pi.ind) > 0){
      
      Pi.draw[Pi.ind] = 1e-10
      
      excess.P = sum(Pi.draw) - 1
      
      ind.max = which.max(Pi.draw)
      Pi.draw[ind.max] = Pi.draw[ind.max] - excess.P
    }
    
    Pi[j, ] = Pi.draw
    Pi.prime.mat[j, ] = Pi.prime[-T0]
  }
  return(list(Pi = Pi, Pi.prime = Pi.prime.mat))
}


# Function to draw samples from the full conditional of Beta
# k.mat : JxT0 matrix containing dish indices for each table in each restaurant
# gam : concentration parameter of the global weights
# L : truncation level for the global weights
update_Beta_2 = function(k.mat, gam, L){
  
  n_prime = sapply(seq_len(L), function(k) sum(k.mat == k))
  n = sum(n_prime)
  
  # update Beta using stick-breaking formulation
  Beta.prime = sapply(seq_len(L-1), function(h)
    rbeta(1, shape1 = (n_prime[h]+1), shape2 = (gam + n - n_prime[h])))
  
  Beta.prime.prod = c(1, cumprod(1- Beta.prime))
  Beta.prime[L] = 1
  
  Beta = sapply(seq_len(L), function(i) Beta.prime[i]*Beta.prime.prod[i])
  
  # setting an lower bound of 10^(-10) for \Beta_k's to avoid numerical issues
  Beta.ind = which(Beta < 1e-10)
  
  if(length(Beta.ind) > 0){
    
    Beta[Beta.ind] = 1e-10
    
    excess.P = sum(Beta) - 1
    
    ind.max = which.max(Beta)
    Beta[ind.max] = Beta[ind.max] - excess.P
  }
  
  return(Beta)
}


# Function to draw samples from the full conditional of alpha0
# Pi.prime : Beta-distributed random variables that make up the stick breaking formulation of Pi.
# (a0, b0) : parameters for the prior distribution of alpha0, alpha0 follows Gamma(a0, b0)
update_alpha0 = function(Pi.prime, a0, b0){
  
  J = nrow(Pi.prime)
  K = ncol(Pi.prime) # K = T0 - 1
  
  alpha0 = rgamma(n = 1, shape = J*K + a0, rate = b0 - sum(log(1-Pi.prime)))
  return(alpha0)
}


# Truncated Unmarginalized Blocked Gibbs Sampler that updates both global and local cluster labels
# x : list of length J, x[[j]] contains data in the jth population
# T.max : truncation level for the local weights
# L.max : truncation level for the global weights
# gam : concentration parameter of global DP
# (a0, b0) : parameters for the prior distribution of alpha0, alpha0 follows Gamma(a0, b0)
# phi.param : parameters specifying the prior distribution of \phi
# Burn.in : Burn in period of the MCMC chain
# M : number of MCMC samples required
# est.density : TRUE/FALSE indicating whether the density needs to be estimated
# y.grid : grid points for estimating the density of the populations, if est.density  = TRUE
u_blocked_gibbs = function(x, T.max, L.max, gam, a0, b0, phi.param,
                           Burn.in, M, est.density = FALSE, y.grid = NULL){
  
  J = length(x)
  
  # set initial values for running the Gibbs sampler
  Pi = matrix(1/T.max, nrow = J, ncol = T.max)
  Beta = rep(1/L.max, L.max)
  
  k.mat = matrix(sample(1:L.max, size = J*T.max, replace = TRUE), nrow = J)
  
  xi = phi.param[1]; lambda = phi.param[2]; tau = phi.param[3]
  phi = rnorm(L.max, mean = xi, sd = 1/sqrt(lambda))
  
  alpha0 = rgamma(1, shape = a0, rate = b0)
  
  # list to store the posterior samples
  Iterates = vector(mode = "list", length = M)
  
  for(m in 1:(M + Burn.in)){
    
    # time at the beginning
    T1 = Sys.time()
    
    # update t
    t.list = update_t(k.mat = k.mat, Pi = Pi, x = x, phi = phi, phi.param = phi.param)
    
    # update k
    k.mat = update_k(t.list = t.list, T0 = T.max, Beta = Beta, x = x, phi = phi, phi.param = phi.param)
    
    # get Z
    Z = get_Z(k = k.mat, t = t.list)
    
    # update phi : same as the function used for BGS
    # BGS.R needs to be sourced for using this function
    phi = update_phi(Z = Z, L = L.max, x = x, phi.param = phi.param)
    
    # update Pi
    res = update_Pi_2(t.list = t.list, alpha0 = alpha0, T0 = T.max)
    Pi = res$Pi
    Pi.prime = res$Pi.prime
    
    # update Beta 
    Beta = update_Beta_2(k.mat = k.mat, gam = gam, L = L.max)
    
    #update alpha0
    alpha0 = update_alpha0(Pi.prime = Pi.prime, a0 = a0, b0 = b0)
    
    # time at the end of all updates
    T2 = Sys.time()
    Tdiff =  difftime(T2, T1, units = "secs")
    
    # ordered phi's
    n.group = sapply(seq_len(L.max), function(j) sum(unlist(Z)==j))
    
    dat = data.frame(phi, n.group, ind = 1:L.max)
    dat = dat[order(-dat$n.group),]
    
    f = NULL
    
    if(est.density == TRUE){
      n.grid = length(y.grid)
      
      # J x n.grid matrix to store the densities along the rows.
      f = matrix(NA, nrow = J, ncol =  n.grid)
      for(j in 1:J){
        
        Pi_j_star = sapply(1:L.max, function(m) sum(Pi[j, (k.mat[j, ] == m)]) )
        
        # evaluate the density for each population
        f[j, ] = sapply(1:n.grid, function(ii) 
          sum( Pi_j_star * dnorm(y.grid[ii], mean = phi, sd = sqrt(1/tau)) ))
        
      }
    }
    
    # print every 200th iteration
    if(m %% 200 == 0){
      print(paste("iteration :", m))
    }
    
    # store samples after Burn in
    if(m > Burn.in){
      Iterates[[m-Burn.in]] = list("Z" = Z, "Pi" = Pi, "phi" = phi, "Beta" = Beta, 
                                   "t" = t.list, "k" = k.mat, "phi.ord" = dat$phi,
                                   "n.group" = dat$n.group, "alpha0" = alpha0,
                                   "indices" = dat$ind, "density" = f, "time" = Tdiff)
    }
  }
  return(Iterates)
  
}