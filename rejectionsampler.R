# Code to sample from the tilted gamma density, 
# f_k(x) \propto (1/\Gamma(x))^J x^{A-1} e^{-B_k x}, x > 0.
# using our proposed rejection sampler

# function to draw n random samples from a truncated gamma distribution
# shape = A, rate = B
# truncation region = (u0, v0)
rtrunc.gamma = function(n, A, B, u0, v0){
  
  F.u0 = pgamma(u0, shape = A, rate = B)
  F.v0 = pgamma(v0, shape = A, rate = B)
  
  # sample using modified inverse cdf technique for truncated distributions
  U = runif(n, min = F.u0, max = F.v0)
  X = qgamma(U, shape = A, rate = B)
  
  return(X)
}

# function to draw n random samples from a truncated exponential distribution
# parameter = lambda
# truncation region = (u0, v0)
rtrunc.exp = function(n, lambda, u0, v0){
  
  # sample using inverse cdf technique
  U = runif(n)
  
  Y = U*exp(lambda*v0) + (1-U)*exp(lambda*u0)
  X = (1/lambda)*log(Y)
  
  return(X)
}

# log density, f_k(.)
# We use independent samples from density, f_k to get samples of beta, top level weight parameters
log.f_k = function(t, J, A, B){
  if(t == Inf){
    res = -Inf
  }
  else res = -(J*lgamma(t)) - (B * t) + ((A-1)*log(t))
  return(res)
}

# derivative of the log density
log.f_k.prime = function(t, J, A, B){
  
  res = -(J*digamma(t)) - B + ((A-1)/t)
  return(res)
}

# Function to find the mode of density f_k
mode.f_k = function(J, A, B, M = 1.5){
  h = function(x){
    return(log.f_k(x, J = J, A = A, B = B))
  }
  m = optimise(h, interval = c(0,M), maximum = TRUE, tol = 0.000001)$maximum
  return(m)
}

# log of un-normalized gamma density for the cover when B>0
log.g_k1 = function(t, J, A, B){
  D1 = digamma(2)
  B1 = (J*D1) + B
  
  res = (2*J*D1) - (t*B1) + ((A-1)* log(t))
  return(res)
}

# log of un-normalized gamma density for the cover when B<0
log.g_k2 = function(t, J, A, B, M){
  gam0 = - digamma(1)
  k1 = (J * gam0*log(M)) - J
  k2 = B + (J*log(M)) - J
  
  res = k1 - (t*k2) + (A-1)*log(t)
  return(res)
}

# equation of tangent line of the log density at point "m"
tangent.eq = function(t, m, a, lambda){
  
  res = a + lambda*(t - m)
  return(res)
}

# Function to draw a sample from the resulting cover density
samp.mixture = function(J, A, B){
  
  # sampler for positive B_k
  if(B > 0){
    
    # get the mode of f_k
    m0 = mode.f_k(J = J, A = A, B = B)
    
    # select the points m1 and m2 at the left and right of the mode
    m1 = m0/2
    m2 = m0 + (1.5 - m0)*0.5
    
    lambda1 = log.f_k.prime(m1, J, A, B)
    lambda2 = log.f_k.prime(m2, J, A, B)
    
    a1 = log.f_k(m1, J, A, B); a2 = log.f_k(m2, J, A, B); a0 = log.f_k(m0, J, A, B)
    
    # points of intersection of the tangent lines
    z1 = (a0 - a1 + m1*lambda1)/lambda1
    z2 = (a0 - a2 + m2*lambda2)/lambda2
    z3 = 1.5
    
    D1 = digamma(2)
    B1 = (J*D1) + B
    c0 = 2*J*D1
    
    # Calculate the ratio of the normalizing constants of the four densities 
    C1byC2 = exp(log(1 - exp(-z1 * lambda1)) - log(lambda1) - log(z2 - z1))
    
    C3byC2 = exp(log( 1 - exp((z3-z2)*lambda2) ) - log(-lambda2) - log(z2 - z1))
    
    C4byC2 = exp(pgamma(z3, shape = A, rate = B1, log.p = TRUE, lower.tail = FALSE) + 
                   lgamma(A) - (A*log(B1)) + c0 - a0 - log(z2 - z1))
    
    C3byC1 = exp(log(lambda1) - log(-lambda2) + log(1 - exp((z3-z2)*lambda2) ) - 
                   log(1 - exp(-z1 * lambda1)))
    
    C4byC1 = exp(pgamma(z3, shape = A, rate = B1, log.p = TRUE, lower.tail = FALSE) +
                   lgamma(A) - (A*log(B1)) + c0 - a0 + log(lambda1) - 
                   log(1 - exp(-z1 * lambda1)))
    
    C4byC3 = exp(pgamma(z3, shape = A, rate = B1, log.p = TRUE, lower.tail = FALSE) +
                   lgamma(A) - (A*log(B1)) + c0 - a0 + log(-lambda2) - 
                   log( 1 - exp((z3-z2)*lambda2) ))
    
    # calculate the weights of the mixture density
    # w1 = C1/C ; w2 = C2/C;  w3 = C3/C
    
    # w1 = C1/(C1+C2+C3+C4) = 1/(1 + C2/C1 + C3/C1 + C4/C1)
    w1 = (1/C1byC2) + C3byC1 + C4byC1
    w1 = 1/(1+w1)
    
    # w2 = C2/(C1+C2+C3+C4) = 1/(1 + C1/C2 + C3/C2 + C4/C2)
    w2 = C1byC2 + C3byC2 + C4byC2
    w2 = 1/(1+w2)
    
    # w3 = C3/(C1+C2+C3+C4) = 1/(1 + C1/C3 + C2/C3 + C4/C3)
    w3 = (1/C3byC1) + (1/C3byC2) + C4byC3
    w3 = 1/(1+w3)
    
    # sample from the mixture density
    U = runif(1)
    
    if(U < w1){
      # sample from g_k1
      S = rtrunc.exp(n = 1, lambda = lambda1, u0 = 0, v0 = z1)
    }
    else if(U < (w1+w2)){
      # sample from g_k2
      S = runif(n = 1, min = z1, max = z2)
    }
    else if(U < (w1+w2+w3)){
      # sample from g_k3
      S = rtrunc.exp(n = 1, lambda = lambda2, u0 = z2, v0 = 1.5)
    }
    else{
      # sample from g_k4
      S = rtrunc.gamma(n = 1, A = A, B = B1, u0 = 1.5, v0 = Inf)
    }
    
    dist = list("M" = 1.5, "m0" = m0, "m1" = m1, "m2" = m2, 
                "z1" = z1, "z2" = z2, "a0" = a0, "a1" = a1, "a2" = a2, 
                "lambda0" = 0, "lambda1" = lambda1, "lambda2" = lambda2)
    
  }
  else{   # sampler for negative B_k
    
    M = exp(((0.01 -B)/J)+1)
    # get the mode of f_k
    m0 = mode.f_k(J = J, A = A, B = B, M = M)
    
    if(m0 < M){
      
      # select the points m1 and m2 at the left and right of the mode
      # if the mode lies in (0, M)
      m1 = m0/2
      m2 = m0 + (M - m0)*0.5
      
      lambda1 = log.f_k.prime(m1, J, A, B)
      lambda2 = log.f_k.prime(m2, J, A, B)
      
      a1 = log.f_k(m1, J, A, B); a2 = log.f_k(m2, J, A, B); a0 = log.f_k(m0, J, A, B)
      
      # points of intersection of the tangent lines
      z1 = (a0 - a1 + m1*lambda1)/lambda1
      z2 = (a0 - a2 + m2*lambda2)/lambda2
      z3 = M
      
      gam0 = - digamma(1)
      B1 = B + (J*log(M)) - J
      c0 = (J * gam0*log(M)) - J
      
      # Calculate the ratio of the normalizing constants of the four densities 
      C1byC2 = exp(log(1 - exp(-z1 * lambda1)) - log(lambda1) - log(z2 - z1))
      
      C3byC2 = exp(log( 1 - exp((z3-z2)*lambda2) ) - log(-lambda2) - log(z2 - z1))
      
      C4byC2 = exp(pgamma(z3, shape = A, rate = B1, log.p = TRUE, lower.tail = FALSE) + 
                     lgamma(A) - (A*log(B1)) + c0 - a0 - log(z2 - z1))
      
      C3byC1 = exp(log(lambda1) - log(-lambda2) + log(1 - exp((z3-z2)*lambda2) ) - 
                     log(1 - exp(-z1 * lambda1)))
      
      C4byC1 = exp(pgamma(z3, shape = A, rate = B1, log.p = TRUE, lower.tail = FALSE) +
                     lgamma(A) - (A*log(B1)) + c0 - a0 + log(lambda1) - 
                     log(1 - exp(-z1 * lambda1)))
      
      C4byC3 = exp(pgamma(z3, shape = A, rate = B1, log.p = TRUE, lower.tail = FALSE) +
                     lgamma(A) - (A*log(B1)) + c0 - a0 + log(-lambda2) - 
                     log( 1 - exp((z3-z2)*lambda2) ))
      
      # calculate the weights of the mixture density
      # w1 = C1/C ; w2 = C2/C;  w3 = C3/C
      
      # w1 = C1/(C1+C2+C3+C4) = 1/(1 + C2/C1 + C3/C1 + C4/C1)
      w1 = (1/C1byC2) + C3byC1 + C4byC1
      w1 = 1/(1+w1)
      
      # w2 = C2/(C1+C2+C3+C4) = 1/(1 + C1/C2 + C3/C2 + C4/C2)
      w2 = C1byC2 + C3byC2 + C4byC2
      w2 = 1/(1+w2)
      
      # w3 = C3/(C1+C2+C3+C4) = 1/(1 + C1/C3 + C2/C3 + C4/C3)
      w3 = (1/C3byC1) + (1/C3byC2) + C4byC3
      w3 = 1/(1+w3)
      
      # sample from the mixture density
      U = runif(1)
      
      if(U < w1){
        # sample from g_k1
        S = rtrunc.exp(n = 1, lambda = lambda1, u0 = 0, v0 = z1)
      }
      else if(U < (w1+w2)){
        # sample from g_k2
        S = runif(n = 1, min = z1, max = z2)
      }
      else if(U < (w1+w2+w3)){
        # sample from g_k3
        S = rtrunc.exp(n = 1, lambda = lambda2, u0 = z2, v0 = M)
      }
      else{
        # sample from g_k4
        S = rtrunc.gamma(n = 1, A = A, B = B1, u0 = M, v0 = Inf)
      }
      
      dist = list("M" = M, "m0" = m0, "m1" = m1, "m2" = m2, 
                  "z1" = z1, "z2" = z2, "a0" = a0, "a1" = a1, "a2" = a2, 
                  "lambda0" = 0, "lambda1" = lambda1, "lambda2" = lambda2)
      
    }
    else{
      
      # if the mode does not lie in (0, M)
      # select three equidistant points in (0,M)
      m1 = M/4; m0 = M/2; m2 = 3*M/4
      
      lambda1 = log.f_k.prime(m1, J, A, B)
      lambda0 = log.f_k.prime(m0, J, A, B)
      lambda2 = log.f_k.prime(m2, J, A, B)
      
      a1 = log.f_k(m1, J, A, B); a2 = log.f_k(m2, J, A, B); a0 = log.f_k(m0, J, A, B)
      
      # points of intersection of the tangent lines
      z1 = (a0 - a1 + m1*lambda1 - m0*lambda0)/(lambda1 - lambda0)
      z2 = (a0 - a2 + m2*lambda2 - m0*lambda0)/(lambda2 - lambda0)
      z3 = M
      
      gam0 = - digamma(1)
      B1 = B + (J*log(M)) - J
      c0 = (J * gam0*log(M)) - J
      
      # Calculate the ratio of the normalizing constants of the four densities 
      C1byC2 = exp(log(1 - exp(-z1 * lambda1)) + log(lambda0) - log(lambda1) - ((z2-z1)*lambda0))
      
      C3byC2 = exp(log( 1 - exp(-(z3-z2)*lambda2) ) - log(1 - exp(-(z2-z1)*lambda0)) +
                     (z3-z2)*lambda2 + log(lambda0) - log(lambda2))
      
      C4byC2 = exp(pgamma(z3, shape = A, rate = B1, log.p = TRUE, lower.tail = FALSE) + 
                     lgamma(A) - (A*log(B1)) + c0 - a0 - (z2 - m0)*lambda0 -
                     log(1- exp(-(z2-z1)*lambda0)))
      
      C3byC1 = exp(log(lambda1) - log(lambda2) + log(1 - exp(-(z3-z2)*lambda2) ) - 
                     log(1 - exp(-z1 * lambda1)) +
                     (z2-z1)*lambda0 + (z3 - z2)*lambda2)
      
      C4byC1 = exp(pgamma(z3, shape = A, rate = B1, log.p = TRUE, lower.tail = FALSE) +
                     lgamma(A) - (A*log(B1)) + c0 - a0 - (z1 - m0)*lambda0 + log(lambda1) - 
                     log(1 - exp(-z1 * lambda1)))
      
      C4byC3 = exp(pgamma(z3, shape = A, rate = B1, log.p = TRUE, lower.tail = FALSE) +
                     lgamma(A) - (A*log(B1)) + c0 - a0 - (z2 - m0)*lambda0 - (z3-z2)*lambda2 +
                     log(lambda2) - log( 1 - exp(-(z3-z2)*lambda2) ))
      
      # calculate the weights of the mixture density
      # w1 = C1/C ; w2 = C2/C;  w3 = C3/C
      
      # w1 = C1/(C1+C2+C3+C4) = 1/(1 + C2/C1 + C3/C1 + C4/C1)
      w1 = (1/C1byC2) + C3byC1 + C4byC1
      w1 = 1/(1+w1)
      
      # w2 = C2/(C1+C2+C3+C4) = 1/(1 + C1/C2 + C3/C2 + C4/C2)
      w2 = C1byC2 + C3byC2 + C4byC2
      w2 = 1/(1+w2)
      
      # w3 = C3/(C1+C2+C3+C4) = 1/(1 + C1/C3 + C2/C3 + C4/C3)
      w3 = (1/C3byC1) + (1/C3byC2) + C4byC3
      w3 = 1/(1+w3)
      
      # sample from the mixture density
      U = runif(1)
      
      if(U < w1){
        # sample from g_k1
        S = rtrunc.exp(n = 1, lambda = lambda1, u0 = 0, v0 = z1)
      }
      else if(U < (w1+w2)){
        # sample from g_k2
        S = rtrunc.exp(n = 1, lambda = lambda0, u0 = z1, v0 = z2)
      }
      else if(U < (w1+w2+w3)){
        # sample from g_k3
        S = rtrunc.exp(n = 1, lambda = lambda2, u0 = z2, v0 = M)
      }
      else{
        # sample from g_k4
        S = rtrunc.gamma(n = 1, A = A, B = B1, u0 = M, v0 = Inf)
      }
      
      dist = list("M" = M, "m0" = m0, "m1" = m1, "m2" = m2, 
                  "z1" = z1, "z2" = z2, "a0" = a0, "a1" = a1, "a2" = a2, 
                  "lambda0" = lambda0, "lambda1" = lambda1, "lambda2" = lambda2)
    }
  }
  
  return(list("sample" = S, "dist" = dist))
}

# log of the cover density
# dist : list containing the parameters characterizing the cover density
log.mixture_dens = function(t, J, A, B, dist){
  
  n = length(t)
  res = rep(0, n)
  
  for(i in 1:n){
    
    if(t[i] == Inf){
      res[i] = -Inf
    }
    else{
      if(t[i] >=0 && t[i] < dist$z1){
        res[i] =  tangent.eq(t = t[i], m = dist$m1, a = dist$a1, lambda = dist$lambda1) 
      }
      else if(t[i] >= dist$z1 && t[i] < dist$z2){
        res[i] =  tangent.eq(t = t[i], m = dist$m0, a = dist$a0, lambda = dist$lambda0) 
      }
      else if(t[i] >= dist$z2 && t[i] < dist$M){
        res[i] =  tangent.eq(t = t[i], m = dist$m2, a = dist$a2, lambda = dist$lambda2) 
      }
      else{
        
        if(B > 0){
          res[i] = log.g_k1(t = t[i], J = J, A = A, B = B) 
        }
        else{
          res[i] = log.g_k2(t = t[i], J = J, A = A, B = B, M = dist$M) 
        }
        
      }
    }
    
    
  }
  return(res)
}
