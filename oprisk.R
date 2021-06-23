# this script puts together many useful functions for operational risk analysis
# ========================================
# section 1: general purpose truncated variable constructor
# ========================================
# pdf for a truncated random variable
dtrunc <- function(x, spec, a, b, log = F, ...) {
	dd <- get(paste("d", spec, sep = ""), mode = "function")
	pd <- get(paste("p", spec, sep = ""), mode = "function")
	if (log == F) {
		return(dd(x, ...)/(pd(b, ...) - pd(a, ...)))
	} else if (log == T) {
		return(dd(x, log = T, ...) - log(pd(b, ...) - pd(a, ...)))
	}
}
# cdf for a truncated random variable
ptrunc <- function(q, spec, a, b, ...) {
	pd <- get(paste("p", spec, sep = ""), mode = "function")
	return((pd(q, ...) - pd(a, ...))/(pd(b, ...) - pd(a, ...)))	
}
# inverse cdf for a truncated random variable
qtrunc <- function(p, spec, a, b, ...) {
	qd <- get(paste("q", spec, sep = ""), mode = "function")
	pd <- get(paste("p", spec, sep = ""), mode = "function")	
	q <- qd(pd(a, ...) + p * (pd(b, ...) - pd(a, ...)), ...)
	return(q)
}
# truncated random variable generator
rtrunc <- function(n, spec, a, b, ...) {
	pd <- get(paste("p", spec, sep = ""), mode = "function")	
	qd <- get(paste("q", spec, sep = ""), mode = "function")	
	u <- runif(n, min = 0, max = 1)
	x <- qd(pd(a, ...) + u * (pd(b, ...) - pd(a, ...)), ...)
	return(x)
}
# ========================================
# section 2: generalized log t and generalized t
# ========================================
# generalized log t
# pdf
dglt <- function(x, location_log, scale_log, df, log = F) {
	tt <- (log(x) - location_log) / scale_log
	if (log ==  F) {
		val <- ifelse(x == 0, 0, dt(tt, df, log = F) * (1/x) * (1/scale_log))
	} else if (log == T) {
		val <- ifelse(x == 0, -Inf, dt(tt, df, log = T) + log(1/x) + log(1/scale_log))	
	}
	return(val)
}
# cdf
pglt <- function(q, location_log, scale_log, df) {
	tt <- (log(q) - location_log) / scale_log
	return(pt(tt, df))
}
# inverse cdf
qglt <- function(p, location_log, scale_log, df) {
	val <- exp(location_log + scale_log * qt(p, df))
	return(val)
}
# random variable generator
rglt <- function(n, location_log, scale_log, df) {
	return(exp(location_log + scale_log * rt(n, df)))
}
# =========================================
# generalized t
# =========================================
# pdf
dgt <- function(x, location, scale, df, log = F) {
	tt <- (x - location) / scale
	if (log ==  F) {
		val <- ifelse(x == 0, 0, dt(tt, df, log = F) * (1/scale))
	} else if (log == T) {
		val <- ifelse(x == 0, -Inf, dt(tt, df, log = T) + log(1/scale))	
	}
	return(val)
}
# cdf
pgt <- function(q, location, scale, df) {
	tt <- (q - location) / scale
	return(pt(tt, df))
}
# inverse cdf
qgt <- function(p, location, scale, df) {
	val <- location + scale * qt(p, df)
	return(val)
}
# random variable generator
rgt <- function(n, location, scale, df) {
	return(location + scale * rt(n, df))
}
# ===========================================
# section 3: asset correlation
# ===========================================
# Estimate Asset Correlation MLE
estimate_ac_ml <- function(pd = NULL, df = list(dd, nn)){
	if(is.null(pd)){
		pd <- df$dd/df$nn
		w <- df$nn/sum(df$nn)
		z <- qnorm(pd)
		r <- wt.var(z, w)/(1 + wt.var(z, w))
		d <- wt.mean(z, w) * sqrt(1 - r)
	} else {
		z <- qnorm(pd)
		r <- var(z)/(1 + var(z))
		d <- mean(z) * sqrt(1 - r)
	}
	return(list(apd = pnorm(d), ac = r))
}
# Estimate Asset Correlation GMM
estimate_ac_gmm <- function(pd = NULL, df = list(dd, nn)){
	if(is.null(pd)){
		pd <- df$dd/df$nn
		w <- df$nn/sum(df$nn)
		p <- wt.mean(pd, w)
		v <- wt.var(pd, w)
	} else {
		p<- mean(pd)
		v <- var(pd)
	}
	d <- qnorm(p)
	e_y2 <- v + p^2
	out <- optimize(f = function(r){
		val <- (e_y2 - pmvnorm(upper = c(d, d), 
		  corr = matrix(c(1, r, r, 1), 2, 2)))^2
		return(val)
	}, interval = c(0, 1))
	return(list(apd = p, ac = out$minimum))
}
# Bootstrap Standard Error
# sim: number of simulations
# N: number of records
# ac: asset correlation estimate
# apd: average PD estimate
# f: method of AC estimation
bootstrap_se_ac <- function(sim, N, ac, apd, f){
	out_ac <- out_apd <- numeric(sim)
	d <- qnorm(apd)
	for(i in seq_len(sim)){
		cond_pd <- pnorm((d - sqrt(ac) * rnorm(N))/(sqrt(1 - ac)))
		out <- f(pd = cond_pd)
		out_ac[i] <- out$ac
		out_apd[i] <- out$apd
	}
	return(list(se_apd = sd(out_apd), se_ac = sd(out_ac)))
}
