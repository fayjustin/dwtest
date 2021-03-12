# This R script for dwtest was modified from the package lmtest.
#
# input weighted model from lm()
# m <- lm (y ~ X, weights = cov)
# cov: sum (# of read counts of both alleles) for intra-specific hybrid strains
# cov: sum (fpkm of two orthologous genes)    for inter-specific hybrid strains
# 
# I will also try normalized coverage, acov, for both intra- and inter- specific hybrid strains.
# 
# usage: 
#     source("weigthed_dwtest.R")
#     m   <- lm(X ~ y, weights = cov)
#     res <- weighted_dw.test(m)
#     autop[i]<-res$p.value
# default alternative: "greater"

# define function
weighted_dw.test <- function(
			formula,
			order.by = NULL,
			alternative = c("greater", "two.sided", "less"),
			iterations  = 15,
			exact       = NULL,
			tol         = 1e-10,
			data        = list(),
			weight_check = NULL
		)
# main
{
	dname       <- paste(deparse(substitute(formula)))
	alternative <- match.arg(alternative)

	# read model, residuals, and weights
	if(!inherits(formula, "formula")) {

		if(!is.null(w <- weights(formula))) {
			weight_check = 1;
		}

		X <- if(is.matrix(formula$x))
			formula$x
		     else model.matrix(terms(formula), model.frame(formula))
		y <- if(is.vector(formula$y))
			formula$y
		     else model.response(model.frame(formula))
		w <- weights(formula)
	}
	else {
		mf <- model.frame(formula, data = data)
		y  <- model.response(mf)
		X  <- model.matrix(formula, data = data)
	}  

	if(!is.null(order.by)) {
		if(inherits(order.by, "formula")) {
			z <- model.matrix(order.by, data = data)
			z <- as.vector(z[,ncol(z)])
		} else {
     			z <- order.by
		}
		X <- as.matrix(X[order(z),])
		y <- y[order(z)]
	}

	n <- nrow(X)
	if(is.null(exact)) exact <- (n < 100)
	k <- ncol(X)
    
	res <- lm.fit(X,y)$residuals

# The original dwtest estimates dw from the residuals, r, as the following:
#         sum(r[i+1] - r[i])^2
#    DW = --------------------
#               sum(r^2)
#
# The weighted dwtest estimates dw from the residuals, r, as the following:
#         sum(r[i+1]w[i+1] - r[i]w[i])^2
#    DW = ------------------------------
#                  sum((wr)^2)
#

	if (is.null(weight_check)) {
		dw <- sum(diff(res)^2)/sum(res^2)         # original
	}
	else {
		dw <- sum(diff(res * w)^2)/sum((res*w)^2) # modified with weight
	}
	Q1 <- chol2inv(qr.R(qr(X)))

	if(n < 3) {
		warning("not enough observations for computing a p value, set to 1")
		pval <- 1
	}
	else {

# There might be a problem when we estimate p-values.
#
#		if (!is.null(weight_check)) {
#			X <- X * w
#		}

		if(exact) {
			A  <- diag(c(1,rep(2, n-2), 1))
			A[abs(row(A)-col(A))==1] <- -1
			MA <- diag(rep(1,n)) - X %*% Q1 %*% t(X)
			MA <- MA %*% A
			ev <- eigen(MA)$values[1:(n-k)]
			if(any(Im(ev)>tol)) warning("imaginary parts of eigenvalues discarded")
			ev <- Re(ev)
			ev <- ev[ev>tol]
			pdw <- function(dw) .Fortran(     # calling fortran from R for faster computation 
				"pan",                    # fortran function, pan, from lmtest package
				as.double(c(dw,ev)),
				as.integer(length(ev)),
				as.double(0),
				as.integer(iterations),
				x=double(1),
				PACKAGE = "lmtest"
				)$x
			pval <- switch(
				alternative,
				"two.sided" = (2*min(pdw(dw), 1-pdw(dw))),
				"less"      = (1 - pdw(dw)),
				"greater"   = pdw(dw)
				)
  
			if(is.na(pval) || ((pval > 1) | (pval < 0))) {
				warning("exact p value cannot be computed (not in [0,1]), approximate p value will be used")
				exact <- FALSE
			}
		}

		if(!exact) {
			if(n < max(5, k)) {
				warning("not enough observations for computing an approximate p value, set to 1")
				pval <- 1        
			}
			else {
				AX <- matrix(as.vector(filter(X, c(-1, 2, -1))), ncol = k)
				AX[1,] <- X[1,] - X[2,]
				AX[n,] <- X[n,] - X[(n-1),]
				XAXQ <- t(X) %*% AX %*% Q1
				P <- 2*(n-1) - sum(diag(XAXQ))
				Q <- 2*(3*n - 4) - 2* sum(diag(crossprod(AX) %*% Q1)) + sum(diag(XAXQ %*% XAXQ))
				dmean <- P/(n-k)
				dvar <- 2/((n-k)*(n-k+2)) * (Q - P*dmean)
				pval <- switch(
					alternative,
					"two.sided" = (2*pnorm(abs(dw-dmean), sd=sqrt(dvar), lower.tail = FALSE)),
					"less"      = pnorm(dw, mean = dmean, sd = sqrt(dvar), lower.tail = FALSE),
					"greater"   = pnorm(dw, mean = dmean, sd = sqrt(dvar))
					)
			}
		}
	}
  
	alternative <- switch(
		alternative,
		"two.sided" = "true autocorrelation is not 0",
		"less"      = "true autocorrelation is less than 0",
		"greater"   = "true autocorrelation is greater than 0"
		)

	names(dw) <- "DW"
	RVAL <- list(
		statistic   = dw,
		method      = "Durbin-Watson autocorrelation test (weighted)",
		alternative = alternative,
		p.value     = pval,
		data.name   = dname
		)
	class(RVAL) <- "htest"
	return(RVAL)
}
