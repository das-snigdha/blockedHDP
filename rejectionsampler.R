# Code to sample from the tilted gamma density, 
# f_k(x) \propto (1/\Gamma(x))^J x^{A-1} e^{-B_k x}, x > 0.
# using our proposed rejection sampler

# function to draw n random samples from a truncated exponential distribution
# parameter = lambda
# truncation region = (u0, v0)
rtrunc.exp = function(n, lambda, u0, v0){
  
  # sample using inverse cdf technique
  U = runif(n)
  
  if(v0 == Inf){
    X = u0 + (1/lambda)*log(1- U)
  }
  else{
    Y = exp(-lambda*(v0 -u0))
    X = v0 + (1/lambda)*log(Y + U*(1-Y))
  }
  
  return(X)
}

# log density, f_k(.)
# We use independent samples from density, f_k to get samples of beta, top level weight parameters
log.f_k = function(t, J, A, B){
  
  res = -(J*lgamma(t)) - (B * t) + ((A-1)*log(t))
  res[which(t==Inf)] = -Inf
  
  return(res)
}

# derivative of the log density
log.f_k.prime = function(t, J, A, B){
  
  res = -(J*digamma(t)) - B + ((A-1)/t)
  return(res)
}

# Function to find the mode of density f_k
mode.f_k = function(J, A, B, M){
  h = function(x){
    return(log.f_k(x, J = J, A = A, B = B))
  }
  m = optimise(h, interval = c(0,M), maximum = TRUE, tol = 0.000001)$maximum
  return(m)
}

# equation of tangent line of the log density at point "m"
tangent.eq = function(t, m, a, lambda){
  
  res = a + lambda*(t - m)
  return(res)
}


# Function to get the points of intersection of the tangent lines at given knot points
# m: vector of knot points
# a : log-density at knot points, m
# lambda : derivative of log-density at knot points, m
intersection.points = function(a, lambda, m){
  K0 = length(m) - 1
  res = sapply(1:K0, function(i) 
    (a[i+1] - a[i] + lambda[i]*m[i] - lambda[i+1]*m[i+1])/(lambda[i] - lambda[i+1]))
  
  return(res)
}

# Function to calculate the ratio of normalizing constants of the densities
# done to stabilize the calculation of the weights of the mixture density, g_k
# C_g_{k,i} / C_g_{k,j} for any supplied i > j
C.ratio = function(i, j, N, J, a, lambda, m, q){
  if(i == (2*N+2)){
    
    s_i = sign(lambda[i]); s_j = sign(lambda[j])
    
    X1 = a[i] - a[j] - (m[i] - q[i])*lambda[i] + (m[j] - q[j])*lambda[j]
    
    X2 = log(s_j*lambda[j]) - log(s_i*lambda[i])
    
    X3 = - log(s_j*(exp((q[j+1] - q[j])*lambda[j]) - 1))
    
    res = exp(X1 + X2 + X3)
  }
  else if(i == (N+1)){
    
    s_j = sign(lambda[j])
    
    X1 = a[i] - a[j] + (m[j] - q[j])*lambda[j] + log(q[i+1] - q[i])
    
    X2 = log(s_j*lambda[j])
    
    X3 = - log(s_j*(exp((q[j+1] - q[j])*lambda[j]) - 1))
    
    res = exp(X1 + X2 + X3)
  }
  else if(j == (N+1)){
    
    s_i = sign(lambda[i])
    
    X1 = a[i] - a[j] - (m[i] - q[i])*lambda[i] - log(q[j+1] - q[j])
    
    X2 = - log(s_i*lambda[i])
    
    X3 = log(s_i*(exp((q[i+1] - q[i])*lambda[i]) - 1))
    
    res = exp(X1 + X2 + X3)
  }
  else {
    
    s_i = sign(lambda[i]); s_j = sign(lambda[j])
    
    X1 = a[i] - a[j] - (m[i] - q[i])*lambda[i] + (m[j] - q[j])*lambda[j]
    
    X2 = log(s_j*lambda[j]) - log(s_i*lambda[i])
    
    X3 = log(s_i*(exp((q[i+1] - q[i])*lambda[i]) - 1)) - 
      log(s_j*(exp((q[j+1] - q[j])*lambda[j]) - 1))
    
    res = exp(X1 + X2 + X3)
  }
  
  if(is.na(res)){
    res = Inf
  }
  return(res)
}

# Function to calculate the weights and parameters of the mixture density
param.mixture = function(J, A, B, N){
  
  K = 2*N + 1
  m = rep(NA, K+1); q = rep(NA, K+1)
  
  if(B > 0){
    # central knot point is the mode
    m[N+1] = mode.f_k(J = J, A = A, B = B, M = 1.5)
    
    # last knot point
    m[K+1] = m[N+1]+1.5
    
  }else {
    M = exp(1 - (B/J))
    # central knot point is the mode
    m[N+1] = mode.f_k(J = J, A = A, B = B, M = M)
    
    # central knot point is the mode
    m[K+1] = M
  }
  
  # first and last knot points
  m[1] = m[N+1]/2; m[K] = (m[N+1]+m[K+1])/2
  
  # (N-1) knot points to the left of the mode
  m[2:N] = seq(m[1], m[N+1], length = (N+1))[2:N]
  
  # (N-1) knot points to the right of the mode
  m[(N+2):(2*N)] = seq(m[N+1], m[K], length = (N+1))[2:N]
  
  # derivative of log-density at vector of knot points, m
  lambda = log.f_k.prime(m, J, A, B)
  
  # log-density at vector of knot points, m
  a = log.f_k(m, J, A, B)
  
  q[1] = 0;
  # points of intersection of the tangent lines
  q[2:(K+1)] = intersection.points(a = a, lambda = lambda, m = m)
  
  
  # get the ratio of the normalizing constants of the piecewise covers
  i_vec = rep((K+1):2, K:1)
  j_vec = unlist(sapply(K:1, function(k) seq(1, k, 1)))
  
  C.ratio.vec = sapply(1:length(i_vec), function(k) 
    C.ratio(i = i_vec[k], j = j_vec[k], N, J, a, lambda, m, q))
  
  C.ratio.list = lapply(1:K, function(k) C.ratio.vec[j_vec == k])
  
  C.ratio.inv.list = lapply(2:(K+1), function(k) 1/C.ratio.vec[i_vec == k])
  
  # weights of the mixture density
  weights = rep(NA, K)
  weights[1] = 1 / (1 + sum(C.ratio.list[[1]]))
  weights[2:K] = sapply(2:K, function(k) 
    1/(1 + sum(C.ratio.list[[k]]) + sum(C.ratio.inv.list[[k-1]] )))
  
  return(list(weights = weights, m = m, a = a, lambda = lambda, q = q))
}


# Function to draw a sample from the resulting mixture density
# weights : weights of the mixture density
# lambda : parameter of the exponential densities in the cover
# q : points of intersection of the tangent lines
samp.mixture = function(weights, lambda, q){
  
  K0 = length(weights)
  cum.wts = cumsum(weights)
  U = runif(1)
  
  pos.U = which(sort(c(cum.wts, U)) == U)
  if(pos.U <= K0){
    S = rtrunc.exp(n = 1, lambda = lambda[pos.U], u0 = q[pos.U], v0 = q[pos.U+1])
  }else 
    S = rtrunc.exp(n = 1, lambda = lambda[pos.U], u0 = q[K0+1], v0 = Inf)
  
  return(S)
  
}

# log of the cover density
# m : vector of knot points
# a : log-density at the knot points
# lambda : parameter of the exponential densities in the cover
# q : points of intersection of the tangent lines
log.mixture_dens = function(t, m, a, lambda, q){
  
  if(t == Inf){
    res = -Inf
  }
  else{
    pos.t = which(sort(c(q, t)) == t)
    res = tangent.eq(t = t, m = m[pos.t-1], a = a[pos.t-1], lambda = lambda[pos.t-1]) 
  }
  
  return(res)
}