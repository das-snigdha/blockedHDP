# Code to implement exact slice sampler (SS) (Amini et al., 2019) for HDP taken from
# https://github.com/aaamini/hdpslicer/blob/master/hdp_inference.R
# with minor modifications (commented with CHANGED/ADDED) for fitting a univariate conjugate Gaussian mixture model.

### HDP inference -- slice sampler ###

# CHANGED: default specification of concentration of global DP beta0 from 3 to 1
# ADDED: prec2_y, prec2_phi and y.grid as inputs of the function
#        prec2_y and prec2_phi denote precision of the likelihood and prior respectively
#        y.grid denotes the grid points for evaluating the density of each population

# hdp_slice_sampler <- function(y, beta0=3, gam0=1, ITRmax=50, Kcap=20, Tcap=20,
#                               doubling_factor=1.5, categorical=F, W=10, cat_prior_alpha=NA, randinit=F) {
hdp_slice_sampler <- function(y, beta0=1, gam0=1, ITRmax=50, Kcap=20, Tcap=20,
                              doubling_factor=1.5, categorical=F, W=10, cat_prior_alpha=NA,
                              randinit=F, prec2_y = 1, prec2_phi = 1, y.grid) {
  # y should be list of length J, where each y[[j]] contains the observations for j^th restaurant, i.e., y[[j]] is a matrix of size  n[j] x d where n[j] is the # of observation or customers in res. j
  # for categorical data W should be given and correct!
  J <- length(y)
  # if (categorical){ # multinomial y[[j]] is n[j] x 1
  #   n <- sapply(y, function(x) length(x))
  # }
  # else { # Gaussian y[[j]] is n[j] x 2 
  #   n <- sapply(y, function(x) dim(x)[1])    
  # }
  if (categorical){ # multinomial y[[j]] is n[j] x 1
    y <- lapply(y, function(x) as.matrix(x))
  }
  
  n <- sapply(y, function(x) dim(x)[1])   
  
  #Tv <- rep(1,J)
  #Kv <- rep(1,J)
  #K <- max(Kv)
  
  #Kcap <- 20  # hard-coded for now
  Tjcap <- rep(Tcap,J) # hard-coded for now
  #Tcap <- Tv
  #Kcap <- K
  if (randinit) {
    tb <- lapply(1:J, function(j) sample(1:Tjcap[j], size=n[j], replace=T))
    kb <- lapply(1:J, function(j) sample(1:Kcap, size=Tjcap[j], replace=T))
    z <-  lapply(1:J, function(j) sapply(1:n[j], function(i) kb[[j]][ tb[[j]][i] ] ))
    
  } else {
    tb <- lapply(1:J, function(j) rep(1,n[j]))
    kb <- lapply(1:J, function(j) rep(1,Tjcap[j]))
    z  <- lapply(n,   function(nj) rep(1,nj))
    
  }
  
  u <- lapply(1:J, function(j) runif(n[j]))
  v <- lapply(1:J, function(j) runif(Tjcap[j]))
  
  kb_old <- kb
  tb_old <- tb
  u_old <- u
  v_old <- v
  z_old <- z
  
  if (categorical) { # multinomial
    require(extraDistr)
    #if (any(is.na(cat_prior_alpha))) cat_prior_alpha <- rep(1/W,W)
    if (any(is.na(cat_prior_alpha))) cat_prior_alpha <- rep(1/W,W)
    
    update_phi <- function(y, z, K) {
      tab <- table( factor(unlist(z), levels=1:K), factor(unlist(y), levels=1:W) )
      cat_post_alpha <- sweep(tab, 2, cat_prior_alpha, '+')
      phi <- rdirichlet(K, cat_post_alpha)
      
      list(phi=phi, cat_post_alpha=cat_post_alpha) 
    }
    
    Ker <- function(y, phi) {
      if (length(y) == 0) return(1)
      phi[y]
    }
    
    
  } else { # Gaussian
    update_phi <- update_phi_gauss
    
    # CHANGED: default value of precision of normal likelihood prec2_y from 5^2 to 1
    # Ker <- function(y, phi, prec2_y=5^2) {
    Ker <- function(y, phi, prec2_y = 1) {
      #require(mvtnorm)
      
      # CHANGED: dmvnorm to dnorm for sampling from a univariate normal distribution
      # dmvnorm(y, mean = phi, sigma=diag(1,2)/prec2_y)
      dnorm(y, mean = phi, sd=sqrt(1/prec2_y))
    }
  }
  
  # z_hist <- list()
  converged <- FALSE  # not used now
  itr <- 1
  
  # ADDED: list to store posterior samples
  Iterates = vector(mode = "list", length = ITRmax)
  
  while (itr <= ITRmax && !converged) {
    
    # ADDED: time at the beginning
    T1 = Sys.time()
    
    # undpate gamma
    Pi = list()     # ADDED: list to store pi_jt using formula \pi_j = sum_{t} \gamma_jt \delta_{k_jt}
    gamp <- list()
    gam <- list()
    T_all <- list()
    Tv <- rep(0,J)
    Tj_overflow <- F
    for (j in 1:J) {
      # update gamma
      g_counts <- beta_counts( tb[[j]], Tjcap[j] )
      gamp[[j]] <- rbeta( Tjcap[j], 1 + g_counts[1,], gam0 + g_counts[2,] )
      gam[[j]] <- stick_break_func(gamp[[j]])
      # plot1(gam[[j]])
      #print(u[[j]])
      #print('-')
      T_all[[j]] <- sapply( 1:n[j], function(i) find_tunc_idx(gam[[j]], u[[j]][i]) )
      Tv[j] <- max(T_all[[j]])
      #print(Tv)
      
      if ( Tv[j] > Tjcap[j] ) {
        Tj_overflow <- T
        Tjcap_old <- Tjcap[j]
        Tjcap[j] <- round( doubling_factor*Tjcap[j] )
        # pad v with enough atoms
        v_old[[j]] <- c( v_old[[j]], runif(Tjcap[j]-Tjcap_old) )
        break
      }
    }
    if (Tj_overflow)  {
      cat('Doubling Tjcap.\n')
      kb <- kb_old
      tb <- tb_old
      u <- u_old
      v <- v_old
      z <- z_old
      next
    }
    # round(gam[[j]],2)
    
    k_counts <- beta_counts( unlist(kb), Kcap )
    betap <- rbeta(Kcap, 1 + k_counts[1,], beta0 + k_counts[2,])
    beta <- stick_break_func( betap )
    K_all <- list()
    Kv <- rep(0,J)
    for (j in 1:J) {
      K_all[[j]] <- sapply( 1:Tjcap[j], function(t) find_tunc_idx(beta, v[[j]][t]) )
      Kv[j] <- max( K_all[[j]] )
    }
    K <- max(Kv)
    if (K > Kcap)  {
      cat('Doubling Kcap.\n')
      Kcap <- round(doubling_factor*Kcap)
      kb <- kb_old
      tb <- tb_old
      u <- u_old
      v <- v_old
      z <- z_old
      next
    }
    
    
    # if we got here, it is safe to save current state as old state
    kb_old <- kb
    tb_old <- tb
    u_old <- u
    v_old <- v
    z_old <- z
    
    # ADDED: pass values of precision of likelihood prec2_y and prior precision prec2_phi
    phi_vec <- update_phi(y, z, Kcap, prec2_y, prec2_phi)$phi
    f_vec <- sapply(1:Kcap, function(k) { function(y) Ker(y, phi_vec[k,]) } )
    
    for (j in 1:J) {      
      
      
      # update k
      for (t in 1:Tjcap[j]) {
        prod_list <- lapply( 1:K_all[[j]][t], function(k) f_vec[[k]]( y[[j]][tb[[j]] == t, ] ) )
        prob_vec <- safe_list_prod(prod_list)
        #prob_vec <- sapply( 1:K_all[[j]][t], function(k) prod( f_vec[[k]]( y[[j]][tb[[j]] == t, ] ) ) )
        kb[[j]][t] <- samplePmf( 1, prob_vec)
      }
      
      # update t
      for (i in 1:n[j]){
        prob_vec <- sapply(1:T_all[[j]][i], function(t) f_vec[[  kb[[j]][t] ]]( y[[j]][i,] ) )
        tb[[j]][i] <- samplePmf( 1, prob_vec)
      }
      
      # update u
      u_upper <- sapply(seq(1,n[j]), function(i) gam[[j]][ tb[[j]][i] ])
      u[[j]] <- runif(n[j], 0, u_upper)
      
      # if (any(is.na(u[[j]]))){
      #   print(">>>>>>")
      #   print('tb[j]')
      #   print(tb[[j]])
      #   print('gam[j]')
      #   print(gam[[j]])
      #   print('u[[j]]')
      #   print( u[[j]] )
      #   print("<<<<<<")
      # }
      
      # update v
      v_upper <- sapply(seq(1,Tjcap[j]), function(t) beta[ kb[[j]][t] ])
      v[[j]] <- runif(Tjcap[j], 0, v_upper)
      
      # update z
      z[[j]] <-  sapply(1:n[j], function(i) kb[[j]][ tb[[j]][i] ] )
      
    } # endfor j
    
    # ADDED: time at the end of all updates
    T2 = Sys.time()
    Tdiff =  difftime(T2, T1, units = "secs")
    
    # ADDED: Calculate \pi_j = sum_{t} \gamma_jt \delta_{k_jt}
    for(j in 1:J){
      Pi[[j]] = sapply(1:Kcap, function(m) sum(gam[[j]][ kb[[j]] == m ]) )
    }
    
    # ADDED: Estimate the density for each population using the posterior samples of phi and Pi
    
    n.grid = length(y.grid)
    
    # J x n.grid matrix to store the densities along the rows.
    f = matrix(NA, nrow = J, ncol =  n.grid)
    for(j in 1:J){
      for(i in 1:n.grid){
        f[j, i] = sum( Pi[[j]] * dnorm(y.grid[i], mean = c(phi_vec), sd = sqrt(1/prec2_y)) )
      }
    }
    
    # ADDED: Order the phi's in decreasing order of cluster occupancy
    n.group = sapply(seq_len(Kcap), function(j) sum(unlist(z)==j))
    dat = data.frame(phi = as.vector(phi_vec), n.group, ind = 1:Kcap)
    dat = dat[order(-dat$n.group),]
    
    # ADDED: store posterior samples 
    Iterates[[itr]] = list("Z" = z, "Pi" = Pi, "phi" = as.vector(phi_vec), "phi.ord" = dat$phi,
                           "n.group" = dat$n.group, "density" = f,  "indices" = dat$ind, 
                           "time" = Tdiff)
    
    itr <- itr + 1
    if (itr %% 200 == 0) {
      
      # CHANGED: Print every 200th iteration and remove printing details of beta
      cat(sprintf("Iteration%6d: ",itr),'\n')
      
      # cat(sprintf("%6d: ",itr),'\n')
      # cat(table(round(beta,2)),"\n")
      # beta_hist <- hist(beta)
      # cat(beta_hist$mids,"\n")
      # cat(beta_hist$counts,"\n")
    }
    
    # z_hist[[itr]] <- z
    
  } # end while
  
  #z
  # z_hist
  
  # CHANGED: Return list of all posterior samples instead of just cluster allocations
  return(Iterates)
} # hdp_slice_sampler


safe_list_prod <- function(prod_list){
  log_sums <- sapply(prod_list, function(x) sum(log(x)))
  sapply(log_sums, function(x) exp(x-max(log_sums)))
}


# CHANGED: Default specification of precision of likelihood and prior, prec2_y and prec2_phi
# update_phi_gauss <- function(y, z, K, prec2_y = 5^2, prec2_phi = 1/(1.5^2)) {
update_phi_gauss <- function(y, z, K, prec2_y = 1, prec2_phi = 1) {
  # require(MASS)
  #require(mvtnorm)
  
  z_flat <- unlist(z)
  y_flat <- do.call(rbind,y)
  z_freq <- tabulate(z_flat,K)
  # sapply(1:K, function (k) sum(z_flat==k))
  
  # CHANGED: Add x1 as y_flat instead of x = (x1, x2) = y_flat to account for univariate samples
  # yz_df <- data.frame(y_flat, z=z_flat) 
  yz_df <- data.frame(x1 = y_flat, z=z_flat) 
  temp <- which(z_freq == 0)
  
  # adding observation (y) zero for missing labels, so that we effectively sample phi from the prior for those with missing labels. 
  
  # CHANGED: x1 = 0 instead of x1=0, x2=0 to account for univariate samples
  # yz_df <- rbind(data.frame(x1=0,x2=0, z=temp), yz_df) # 
  yz_df <- rbind(data.frame(x1=0, z=temp), yz_df) # 
  
  phi_mean <- aggregate(.~z, yz_df, function(x) sum(x)*(prec2_y/(prec2_phi + length(x)*prec2_y)) )  # mean of phi-posterior 
  phi_mean <- phi_mean[,-1]
  prec2_post <- prec2_phi + z_freq*prec2_y
  
  # CHANGED: mvrnorm to rnorm to sample from univariate normal distribution
  # phi_new <- do.call(rbind, lapply(1:K, function(k) mvrnorm(1, mu = phi_mean[k,], Sigma = diag(1,2)/prec2_post[k]) ) )
  phi_new <- do.call(rbind, lapply(1:K, function(k) 
    rnorm(1, mean = phi_mean[k], sd = sqrt(1/prec2_post[k])) ) )
  
  list(phi=phi_new, prec2=prec2_post)
}

find_tunc_idx <- function(beta, threshold) {
  # which(cumsum(beta) > 1-threshold)[1]
  # we return the index after the length of beta if there is insufficient atoms
  temp <- which(c(cumsum(beta),1) > 1-threshold)[1]
  if (any(is.na(temp))){
    print('---')
    print(beta)
    print(threshold)
  }
  
  temp
  
}

### Auxiliary functions ###
# samplePmf <- function(n, pmf, normalize=F){
#   if (normalize) { 
#     pmf <- pmf / sum(pmf) 
#   }
#   sample(x = seq(1,length(pmf)), n, replace = T, prob=pmf)  
# }
samplePmf <- function(n, pmf){
  # automatic normalization of pmf by "sample" (?)
  # samples from a single pmf
  sample(x = seq(1,length(pmf)), n, replace = T, prob=pmf)  
}


stick_break_func <- function(x) {
  temp <- c(1,cumprod(1-x))
  
  temp[1:length(x)] * x
}

beta_counts <- function(z,K) {
  # z is a vector of labels, 
  # K is the maximum label number to be counted in z
  # output[1,j]: hom many labels == j are in z
  # output[2,j]: hom many labels > j are in z
  zcounts <- tabulate(z,K)
  
  rbind(zcounts, c(rev(cumsum(rev(zcounts)))[-1],0))    
}