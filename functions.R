library(zoo)
library(sandwich)

############# Omega function #############

meatSQ = meatHC
cpLine = which(grepl("crossprod\\(rval\\)", body(meatSQ)))
# order might be all messed up
body(meatSQ)[[cpLine]] = substitute(rval)
sqrtLine = which(grepl("rval <- sqrt\\(omega\\) \\* X", body(meatSQ)))
body(meatSQ)[[sqrtLine]] = substitute(rval <- list(omega, X))

switchLine = which(grepl("switch", body(meatSQ)))
body(meatSQ)[[switchLine]] = substitute({if (is.null(omega)) {
  type <- match.arg(type)
  if (type == "HC") 
    type <- "HC0"
  switch(type, const = {
    omega <- NULL
  }, HC0 = {
    omega <- function(residuals, diaghat, df) residuals
  }, HC1 = {
    omega <- function(residuals, diaghat, df) residuals * sqrt(length(residuals)/df)
  }, HC2 = {
    omega <- function(residuals, diaghat, df) residuals/sqrt(1 - diaghat)
  }, HC3 = {
    omega <- function(residuals, diaghat, df) residuals/(1 - diaghat)
  }, HC4 = {
    omega <- function(residuals, diaghat, df) {
      n <- length(residuals)
      p <- as.integer(round(sum(diaghat), digits = 0))
      delta <- pmin(4, n * diaghat/p)
      residuals/(1 - diaghat)^(delta/2)
    }
  }, HC4m = {
    omega <- function(residuals, diaghat, df) {
      gamma <- c(1, 1.5)
      n <- length(residuals)
      p <- as.integer(round(sum(diaghat), digits = 0))
      delta <- pmin(gamma[1], n * diaghat/p) + pmin(gamma[2], n * diaghat/p)
      residuals/(1 - diaghat)^(delta/2)
    }
  }, HC5 = {
    omega <- function(residuals, diaghat, df) {
      k <- 0.7
      n <- length(residuals)
      p <- as.integer(round(sum(diaghat), digits = 0))
      delta <- pmin(n * diaghat/p, pmax(4, n * k * 
                                          max(diaghat)/p))
      residuals/((1 - diaghat)^delta)^0.25
    }
  } )
  if (type %in% c("HC2", "HC3", "HC4", "HC4m", "HC5")) {
    if (inherits(diaghat, "try-error")) stop(sprintf("hatvalues() could not be extracted successfully but are needed for %s", 
                                                     type))
    id <- which(diaghat > 1 - sqrt(.Machine$double.eps))
    if (length(id) > 0L) {
      id <- if (is.null(rownames(X))) 
        as.character(id)
      else rownames(X)[id]
      if (length(id) > 10L) 
        id <- c(id[1L:10L], "...")
      warning(sprintf("%s covariances become numerically unstable if hat values are close to 1 as for observations %s", 
                      type, paste(id, collapse = ", ")))
    }
  }
} })


############# Weight function ###############

# Residual multipliers for HC* types (not squared)
# sqrt(omega_i) = |res_i| * g_i  for HC0..HC5
# For "const", sqrt(omega_i) is a constant row scale (not multiplying residuals).
#
# Return:
#   list(weights = numeric(n), kind="residual_multiplier"/"constant_scale", details=list(...))
#
get_residual_multiplier <- function(model,
                                    type = c("HC3","const","HC","HC0","HC1","HC2","HC4","HC4m","HC5"),
                                    treat_const_as_one = FALSE) {
  type <- match.arg(type)
  if (type == "HC") type <- "HC0"
  
  # X, n, p, df with alias protection (like sandwich)
  X <- model.matrix(model)
  if (any(alias <- is.na(coef(model)))) X <- X[, !alias, drop = FALSE]
  n <- NROW(X); p <- NCOL(X); df <- n - p
  
  # hat values + numeric guard
  h <- try(stats::hatvalues(model), silent = TRUE)
  if (inherits(h, "try-error") || anyNA(h)) {
    stop(sprintf("hatvalues() could not be extracted successfully but are needed for %s.", type))
  }
  eps <- .Machine$double.eps
  one_m_h <- pmax(1 - h, eps)
  
  # p estimate used in HC4/HC4m/HC5 (trace of hat matrix rounded, like sandwich)
  p_int <- max(1L, as.integer(round(sum(h), digits = 0)))
  n_over_p <- n / p_int
  out_kind <- "residual_multiplier"
  delta <- NULL
  
  # build multipliers
  if (type == "HC0") {
    w <- rep(1, n)
    
  } else if (type == "HC1") {
    w <- rep(sqrt(n / max(df, 1L)), n)  # guard df<=0
    
  } else if (type == "HC2") {
    w <- 1 / sqrt(one_m_h)
    
  } else if (type == "HC3") {
    w <- 1 / one_m_h
    
  } else if (type == "HC4") {
    # delta_i = pmin(4, n*h_i/p)
    delta <- pmin(4, n_over_p * h)
    w <- 1 / (one_m_h^(delta / 2))
    
  } else if (type == "HC4m") {
    # gamma = (1, 1.5); delta_i = pmin(1, n*h_i/p) + pmin(1.5, n*h_i/p)
    nhp <- n_over_p * h
    delta <- pmin(1, nhp) + pmin(1.5, nhp)
    w <- 1 / (one_m_h^(delta / 2))
    
  } else if (type == "HC5") {
    # k=0.7; delta_i = pmin(n*h_i/p, pmax(4, n*k*max(h)/p))
    k <- 0.7
    delta_cap <- max(4, n_over_p * k * max(h))
    delta <- pmin(n_over_p * h, delta_cap)
    # omega_i = e_i^2 / sqrt((1-h_i)^{delta_i})
    # sqrt(omega_i) = |e_i| / (1 - h_i)^{delta_i/4}
    w <- 1 / (one_m_h^(delta / 4))
    
  } else if (type == "const") {
    if (treat_const_as_one) {
      # 如果只想返回 1（不改变残差），设置 treat_const_as_one=TRUE
      w <- rep(1, n)
      out_kind <- "constant_scale"
    } else {
      # 与 sandwich 保持一致：sqrt(omega) 为常数缩放
      # GLM: working residuals * working weights，并对非 Poisson/Binomial/NB 做 rescale
      if (inherits(model, "glm")) {
        res <- as.vector(residuals(model, type = "working")) * stats::weights(model, type = "working")
        fam <- substr(model$family$family, 1L, 17L)
        if (!(fam %in% c("poisson", "binomial", "Negative Binomial"))) {
          res <- res * sum(stats::weights(model, "working"), na.rm = TRUE) / sum(res^2, na.rm = TRUE)
        }
      } else if (inherits(model, "lm")) {
        res <- as.vector(residuals(model))
        if (!is.null(stats::weights(model))) res <- res * stats::weights(model)
      } else {
        res <- as.vector(residuals(model))
      }
      const_scale <- sqrt(sum(res^2, na.rm = TRUE) / max(df, 1L))
      w <- rep(const_scale, n)
      out_kind <- "constant_scale"
    }
  } else {
    stop("Unhandled type.")
  }
  
  list(
    weights = as.numeric(w),
    kind    = out_kind,
    details = list(
      n = n, p = p, df = df,
      hatvalues = as.numeric(h),
      delta = if (!is.null(delta)) as.numeric(delta) else NULL
    )
  )
}

get_residual_multiplier_vec <- function(model,
                                        type = c("HC3","const","HC","HC0","HC1","HC2","HC4","HC4m","HC5"),
                                        treat_const_as_one = FALSE) {
  get_residual_multiplier(model, type, treat_const_as_one)$weights
}



############# Other ############

resiParam = function(mod, alpha=0.05, type = type, ...){
  vars = labels(terms(formula(mod)))
  # Need to add in multiplying RESI by DoF
  # anova(mod)$coef$Df
  rval = as.data.frame(do.call(rbind, lapply(vars, function(var){ 
    resiPEse(mod, variable=var, alpha=alpha, type = type,...=...)})))
  rownames(rval) = vars
  rval
}


var_sqrtn_beta_hat <- function(alpha, beta, px) {
  pi0 <- exp(alpha) / (1 + exp(alpha))
  pi1 <- exp(alpha + beta) / (1 + exp(alpha + beta))
  
  # var(sqrt(n) betahat)
  S0 <- px * pi1 * (1 - pi1) + (1 - px) * pi0 * (1 - pi0)
  S1 <- px * pi1 * (1 - pi1)
  variance <- S0 / (S0 * S1 - S1^2)
  
  return(variance)
}

# Define S^2 function
S_squared <- function(beta, alpha, px) {
  beta^2 / var_sqrtn_beta_hat(alpha, beta, px)
}


vect <- function(mat){
  matrix(matrixcalc::vec(as.matrix(mat)), nrow = 1)
}

sym <- function(M) (M + t(M))/2


d_to_resi <- function(d, a, beta, px = 0.5, signed = TRUE) {
  # d: Cohen's d (obtained from Chinn's logOR -> d conversion)
  # a: intercept from the logistic model
  # beta: slope coefficient from the logistic model
  # px: probability that X = 1 (default = 0.5)
  # signed: whether to preserve the sign of d (TRUE) or return absolute value (FALSE)
  
  # Step 1: compute probabilities at x=0 and x=1
  mu0 <- plogis(a)          # P(Y=1 | X=0)
  mu1 <- plogis(a + beta)   # P(Y=1 | X=1)
  
  # Step 2: compute Bernoulli variances
  w0  <- mu0 * (1 - mu0)
  w1  <- mu1 * (1 - mu1)
  
  # Step 3: compute the constant scaling factor
  # This accounts for the intercept–slope correlation in the Fisher information
  k <- sqrt(px * (1 - px) * w0 * w1 / ((1 - px) * w0 + px * w1)) / (sqrt(3) / pi)
  
  # Step 4: return RESI equivalent of d
  if (!signed) {
    return(abs(d) * k)
  } else {
    return(d * k)
  }
}

var_Shat_null <- function(m1) {
  if (!requireNamespace("gsl", quietly = TRUE)) {
    stop("Please install the 'gsl' package via install.packages('gsl').")
  }
  
  gamma <- base::gamma
  
  # ---- (1) E[ sqrt{n} S_hat ]
  E1 <- (gamma(3/2) / (2^(m1/2) * gamma(m1/2))) *
    m1^((m1 + 1)/2) *
    gsl::hyperg_U(3/2, -(m1 + 1)/2, m1/2)
  
  # ---- (2) E[ n S_hat^2 ] = E[(chi^2 - m1)_+]
  # gsl::gamma_inc(a, x) gives the *lower* incomplete gamma γ(a, x)
  # so upper incomplete gamma = Γ(a) - γ(a, x)
  Gamma_full <- gamma(m1/2)
  gamma_lower1 <- gsl::gamma_inc(m1/2 + 1, m1/2)
  gamma_lower2 <- gsl::gamma_inc(m1/2, m1/2)
  Gamma_upper1 <- Gamma_full * (1 - gamma_lower1 / Gamma_full)
  Gamma_upper2 <- Gamma_full * (1 - gamma_lower2 / Gamma_full)
  
  E2 <- (2 * Gamma_upper1 / Gamma_full) -
    (m1 * Gamma_upper2 / Gamma_full)
  
  # ---- (3) Variance
  var_val <- E2 - E1^2
  return(var_val)
}

########## CI construction ##########

truncate_ci <- function(Shat, se, n, m1, alpha  = 0.05){
  
  z1 <- qnorm(1 - alpha/2)

  # Normal lower bound
  SL <- Shat - z1 * se
  # Normal upper bound
  SU <- Shat + z1 * se
  
  # --- decide regime ---
  if (SL > 0) { # positive lower limit
    return(c(SL, SU))
    
  }else{
    
    # Calculate the percentile of shat under the null
    beta <- pchisq(m1 + n * Shat^2, df = m1, lower.tail = FALSE)
    
    if(beta < alpha/2){
      # Calculate the new lower bound
      SU <- Shat + qnorm(1- (alpha - beta)) * se
      
      return(c(0, SU))
      
    }else{
      SU <- Shat + qnorm(1 - alpha) * se
      
      return(c(0, SU))
    }
  }
}



########## ANOVA ##########

# get_Anova_L():  Return the linear contrast matrix L
#                 used by car::Anova(type = 2 or 3) for a given term

get_L_anova2 <- function(mod, term, type = c("HC3", "const", "HC", "HC0", "HC1", "HC2", 
                                             "HC4", "HC4m", "HC5"), ...) {
  
  # ---------- Extract term information ----------
  type <- match.arg(type)
  names <- labels(terms(mod))
  which.term <- which(term == names)
  not.aliased <- !is.na(coef(mod))
  factors <- attr(terms(mod), "factors")
  
  # assign vector (works for both lm and glm)
  assign <- attr(model.matrix(mod), "assign")
  assign[!not.aliased] <- NA
  subs.term <- which(assign == which.term)
  
  # ---------- Skip if single term ----------
  if (length(names) == 1) 
    return(NULL)
  
  # ---------- Identify relatives ----------
  relatives <- (1:length(names))[-which.term][
    sapply(names[-which.term], function(term2) all(factors[, term] <= factors[, term2]))]
  
  subs.relatives <- unlist(lapply(relatives, function(rel) {
    which(assign == rel)
  }))

  # ---------- Identity matrix ----------
  I.p <- diag(length(coefficients(mod)))
  
  hyp.matrix.1 <- I.p[subs.relatives, , drop = FALSE]
  hyp.matrix.1 <- hyp.matrix.1[, not.aliased, drop = FALSE]
  hyp.matrix.2 <- I.p[c(subs.relatives, subs.term), , drop = FALSE]
  hyp.matrix.2 <- hyp.matrix.2[, not.aliased, drop = FALSE]
  
  # ---------- Use custom covariance ----------
  if(inherits(mod, "geeglm")){
    Sigma = mod$geese$vbeta * n_distinct(mod$id)
  }else{
    if(type == "const"){
      Sigma <- vcov(mod, ...)
    }else{
      Sigma <- vcovHC(mod, type, ...)
    }
  }
  

  if (nrow(hyp.matrix.1) == 0){
    hyp.matrix.term <- hyp.matrix.2
  }else{
    hyp.matrix.term <- t(conjcomp(
      t(hyp.matrix.1),
      t(hyp.matrix.2),
      Sigma
    ))
  }
  
  # Remove all-zero rows
  hyp.matrix.term <- hyp.matrix.term[!apply(hyp.matrix.term, 
                                            1, function(x) all(x == 0)), , drop = FALSE]

  return(hyp.matrix.term)
}


conjcomp <- function (X, Z = diag(nrow(X)), ip = diag(nrow(X))) 
{
  xq <- qr(t(Z) %*% ip %*% X)
  if (xq$rank == 0) 
    return(Z)
  Z %*% qr.Q(xq, complete = TRUE)[, -(1:xq$rank)]
}


########## Final function ##########

resiPEse = function(mod, alpha = 0.05, type = type){
  if(inherits(mod, "geeglm")){
    type = "HC0"
  }
  type <- match.arg(type)
  
  # ---------- 1. ANOVA-level RESI ----------
  vars_anova <- labels(terms(mod))
  resi_anova <- do.call(rbind, lapply(vars_anova, function(term) {
    L <- get_L_anova2(mod, term, type = type)
    if(inherits(mod, "geeglm")){
      resi_a_geeglm(mod, variable = term, alpha = alpha, 
                    unsigned = TRUE)
    }else{
      resiPEse_Th1(mod, L = L, variable = term, alpha = alpha, type = type, 
                   unsigned = TRUE)
    }
  }))
  rownames(resi_anova) <- vars_anova
  
  # ---------- 2. Coefficient-level RESI ----------
  vars_coef <- names(coef(mod))
  resi_coef <- do.call(rbind, lapply(vars_coef, function(var) {
    if(inherits(mod, "geeglm")){
      resi_a_geeglm(mod, variable = var, alpha = alpha, unsigned = FALSE)
    }else{
      resiPEse_Th1(mod, variable = var, alpha = alpha, type = type, unsigned = FALSE)
    }
  }))
  rownames(resi_coef) <- vars_coef
  
  list(anova = resi_anova, coefficients = resi_coef)
}
