resiPEse_Th2 <- function(model, variable = NULL, L = NULL,
                      type = c("HC3", "const", "HC", "HC0", "HC1", "HC2", "HC4", "HC4m", "HC5"),
                      alpha = 0.05, unsigned = FALSE) {
  
  type <- match.arg(type)
  X_full <- model.matrix(model)
  if (any(alias <- is.na(coef(model)))) X_full <- X_full[, !alias, drop = FALSE]
  
  ef <- sandwich::estfun(model)
  ef <- ef[, colnames(X_full), drop = FALSE]
  
  if(!is.null(variable)){
    if(unsigned){
      idx <- if(grepl(":", variable)) which(startsWith(colnames(X_full), variable)) else 
        which(colnames(X_full) == variable | (startsWith(colnames(X_full), variable) & !grepl(":", colnames(X_full))))
    }else{
      idx <- which(colnames(X_full) == variable)
    }
    m1 <- length(idx)
  }
  
  p <- ncol(X_full)
  n <- nrow(X_full)
  
  e <- residuals(model, "response")
  X <- X_full
  
  ## ========== Linear model ==========
  
  if (class(model)[1] == "lm") {
    
    # Add dispersion phi as first element of theta
    phi <- summary(model)$sigma^2
    m <- p + 1
    
    # Estimating equations with phi included
    psi_list <- lapply(1:n, function(i) {
      xi <- X[i, , drop = FALSE]
      ei <- e[i]
      
      psi_phi <- (ei^2 - phi) / (2 * phi^2)
      psi_beta <- (ei * xi) / phi
      cbind(psi_phi, psi_beta)
    })
    
    # First derivative of psi
    psiprime_list <- lapply(1:n, function(i) {
      xi <- X[i, , drop = FALSE]   # 1 x p
      ei <- e[i]
      
      d_phi_phi <- (phi - 2 * ei^2) / (2 * phi^3)  # scalar
      d_phi_beta <- -ei * xi / phi^2                # 1 x p
      d_beta_phi <- t(d_phi_beta)                   # p x 1
      d_beta_beta <- -crossprod(xi) / phi           # p x p
      
      rbind(c(d_phi_phi, d_phi_beta),               # 1 x (1+p)
            cbind(d_beta_phi, d_beta_beta))         # p x (1+p)
    })
    
    # ========== GLM ==========
    
  } else if (class(model)[1] == "glm") {
    
    m <- p
    
    mu_hat <- pmin(pmax(fitted(model), .Machine$double.eps), 1 - .Machine$double.eps)
    w_vec <- weights(model, type = "working")
    # psi_list <- lapply(1:n, function(i) e[i] * X[i, , drop = FALSE]) # Psi_i
    psi_list <- lapply(1:n, function(i) matrix(ef[i, , drop = TRUE], nrow = 1))
    psiprime_list <- lapply(1:n, function(i) -w_vec[i] * crossprod(X[i, , drop = FALSE]))
    
    if(model$family$family %in% c("binomial", "poisson")){
      phi <- 1
    }else{
      phi <- sum((model$y - mu_hat)^2 / w_vec) / model$df.residual
    }
    
    
  } else stop("Unsupported model type")
  
  # ========== Coefficients ==========
  
  if(is.null(L)){
    L <- matrix(0, nrow = m1, ncol = m)
    L_model <- matrix(0, nrow = m1, ncol = p)
    
    if(class(model)[1] == "lm"){
      for (i in 1:m1) {
        L[i, 1 + idx[i]] <- 1
      }
      for (i in 1:m1) {
        L_model[i, idx[i]] <- 1
      }
    }else{
      for (i in 1:m1) {
        L[i, idx[i]] <- 1
      }
      L_model <- L
    }
  }else{
    if(class(model)[1] == "lm"){
      L_model <- L
      L <- as.matrix(cbind(0, L_model))
    }else{
      L_model <- L
    }
  }
  
  if(class(model)[1] == "lm"){
    theta <- c(phi, coef(model))
    
  }else{
    theta <- coef(model)
  }
  
  m1 <- nrow(L)
  beta <- L %*% theta
  
  # ========== Taylor Expansion ==========
  
  # All in matrix form 
  A_full <- - Reduce("+", psiprime_list)/n 
  A_full <- sym(A_full)
  A_inv <- tryCatch(
    chol2inv(chol(A_full)),
    error = function(e) solve(A_full)
  )

  
  ## Matrix B
  h        <- hatvalues(model)
  one_m_h  <- pmax(1 - h, .Machine$double.eps)
  hc_sqrtw <- switch(toupper(type),
                     "HC0"  = rep(1, n),
                     "HC1"  = rep(sqrt(n / max(n - ncol(X), 1L)), n),
                     "HC2"  = 1 / sqrt(one_m_h),
                     "HC3"  = 1 / (one_m_h),
                     "HC4"  = {
                       p_int <- max(1L, as.integer(round(sum(h), 0)))
                       nhp   <- (n / p_int) * h
                       one_m_h^(-pmin(4, nhp) / 2)
                     },
                     "HC4M" = {
                       p_int <- max(1L, as.integer(round(sum(h), 0)))
                       nhp   <- (n / p_int) * h
                       delta <- pmin(1, nhp) + pmin(1.5, nhp)
                       one_m_h^(-delta / 2)
                     },
                     "HC5"  = {
                       p_int    <- max(1L, as.integer(round(sum(h), 0)))
                       nhp      <- (n / p_int) * h
                       k        <- 0.7
                       deltaCap <- max(4, (n / p_int) * k * max(h))
                       delta    <- pmin(nhp, deltaCap)
                       one_m_h^(-delta / 4)
                     },
                     rep(1, n)
  )
  sqrtw <- if (type == "const") rep(1, n) else hc_sqrtw
  cfac  <- if (toupper(type) == "HC1") n/(n-p) else 1
  
  B_full <- cfac * Reduce("+", lapply(1:n, function(i) {
    psi_i <- psi_list[[i]]        # (1 x m)
    if (class(model)[1] == "lm") {
      psi_i[1, 2:m] <- psi_i[1, 2:m] * sqrtw[i]    # lm: only beta 
    }else{
      psi_i[1, 1:m] <- psi_i[1, 1:m] * sqrtw[i]
    }
    crossprod(psi_i)                                # (m x m)
  })) / n
  B_full <- sym(B_full)
  
  # Covariance
  if("glm" %in% class(model) & type == "const"){
    cov_theta_full <- sym(A_inv)
    model_cov_full <- sym(vcov(model)*n)
  }else{
    cov_theta_full <- sym(A_inv %*% B_full %*% A_inv)
    model_cov_full <- sym(vcovHC(model, type)*n)
  }
  
  # All in scalar 
  cov_beta <- L %*% cov_theta_full %*% t(L)
  model_cov_beta <- L_model %*% model_cov_full %*% t(L_model)
  if(m1>1){
    diag(model_cov_beta) <- pmax(diag(model_cov_beta), 0)
  }
  
  
  # ========== Derivatives ==========
  
  
  # # --- theta term --- 
  theta_diff <- lapply(psi_list, function(i) A_inv %*% t(i))
  # # --- A term --- 
  A_diff_vec <- lapply(1:n, function(i) vec(- psiprime_list[[i]] - A_full))
  # # --- B term --- 
  B_diff_vec <- lapply(1:n, function(i) {
    psi_i <- psi_list[[i]]        # (1 x m)
    if (class(model)[1] == "lm") {
      psi_i[1, 2:m] <- psi_i[1, 2:m] * sqrtw[i]    
    } else { # glm
      psi_i[1, ] <- psi_i[1, ] * sqrtw[i]   
    }
    vec(crossprod(psi_i) - B_full)                                # (m x m)
  })
  
  deriv_theta <- t(beta) %*% solve(cov_beta) %*% L
  
  if(type == "const"){
    deriv_A <-  vect(tcrossprod(solve(cov_beta) %*% beta)) %*%
      kronecker(L %*% A_inv, L%*% A_inv) / 2
  }else{
    deriv_A <- vect(tcrossprod(solve(cov_beta) %*% beta)) %*% 
      (kronecker(L %*% A_inv, L %*% cov_theta_full) +
         kronecker(L %*% cov_theta_full, L %*% A_inv)) / 2
  }
  
  deriv_B <- - vect(tcrossprod(solve(cov_beta) %*% beta)) %*% 
    kronecker(L %*% A_inv, L %*% A_inv) / 2
  

  term_theta <- sapply(theta_diff, function(i) deriv_theta %*% i) / sqrt(n)
  term_A <- sapply(A_diff_vec, function(i) deriv_A %*% i)/ sqrt(n)
  term_B <- sapply(B_diff_vec, function(i) deriv_B %*% i) / sqrt(n)
  
  
  if(type == "const"){
    psis <- rbind(term_theta, term_A)
  }else{
    psis <- rbind(term_theta, term_A, term_B)
  }
  
  
  # ========== Final RESI with SE ========== 
  
  vars <- tcrossprod(psis)
  Ssq <- t(beta) %*% solve(model_cov_beta) %*% beta
  
  quad <- sum(vars)
  quad <- pmax(quad, 0)
  
  Z <- beta %*% (model_cov_beta/n)^(-1/2)
  resiSE <- if (Ssq <= 0) 1 / sqrt(n) else sqrt(quad / Ssq / n)
  
  # Xinyu RESI
  if(unsigned){
    if(class(model)[1] == "lm"){
      resi_xinyu <- f2S(n * Ssq/m1, df = m1, rdf = n - p, n = n)
    }else{
      resi_xinyu <- chisq2S(n * Ssq, df = m1, n = n)
    }
    
  }else{
    if(class(model)[1] == "lm"){
      resi_xinyu <- t2S(Z, rdf = model$df.residual, n = n, unbiased = TRUE)
    }else{
      resi_xinyu <- z2S(Z, n = n, unbiased = TRUE)
    }
  }
  
  # Megan RESI
  if("glm" %in% class(model) & type == "const"){
    if(unsigned){
      resi_megan <- RESI::resi_pe(model, vcovfunc = vcov, 
                                  unbiased = TRUE)$anova[variable, "RESI"]
    }else{
      resi_megan <- RESI::resi_pe(model, vcovfunc = vcov, 
                                  unbiased = TRUE)$coefficients[variable, "RESI"]
    }
    
  }else{
    if(unsigned){
      resi_megan <- RESI::resi_pe(model, vcov.args = list(type = type), 
                                  unbiased = TRUE)$anova[variable, "RESI"]
    }else{
      resi_megan <- RESI::resi_pe(model, vcov.args = list(type = type), 
                                  unbiased = TRUE)$coefficients[variable, "RESI"]
    }
  }
  
  # Results
  if(unsigned){
    rval <- data.frame(Df = m1, RESI_Xinyu = resi_xinyu, RESI_Megan = resi_megan, RESI.SE = resiSE)
    rval[, c("LCI", "UCI")] <- truncate_ci(Shat = resi_xinyu, se = resiSE, n = n, m1 = m1)
    
  }else{
    rval <- data.frame(Estimate = beta, Std.Error = sqrt(cov_beta/n),
                       RESI_Xinyu = resi_xinyu, RESI_Megan = resi_megan, RESI.SE = resiSE)
    rval[, c("LCI", "UCI")] <-
      cbind(
        as.vector(resi_xinyu) - qnorm(alpha / 2, lower.tail = FALSE) * rval[, "RESI.SE"],
        as.vector(resi_xinyu) + qnorm(alpha / 2, lower.tail = FALSE) * rval[, "RESI.SE"]
      )
  }
  
  rval
}