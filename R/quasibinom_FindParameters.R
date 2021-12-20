#' Find Quasi-Binomial Distribution Parameters from Moments
#'
#' @description Find the parameters `prob` and `phi` for the quasibinomial
#' given a desired expectation, variance, and size.
#'
#' @param size number of trials (greater than 0).
#' @param mean desired mean of a quasibinomial distribution.
#' @param var desired variance of a quasibinomial distribution.
#'
#' @return A numeric vector. `prob` is probability of success on each trial only if phi = 0.
#' `phi` is the dispersion parameter. Distribution is equivalent to a binomial when phi is 0.
#'
#' @export
#'
#' @details
#' p(x) =  	\choose{\text{size},k} \cdot \text{prob} \cdot
#'          (\text{prob} + x\cdot\phi)^(x-1) \cdot (1-\text{prob}-x\cdot\phi)^(\text{size}-x)
#'
#' Phi must be <= min(prob/size, 1 - prob/size) for a proper density function.
#' The function `quasibinom_FindParameters` uses optim to find the parameter
#' combination to match a given expectation and variance. Answer should not be
#' interpreted as exact, so proceed with caution.
#'
#' @examples
#'
#' Equasibinom(30, 0.5, 0.1/30)
#' Vquasibinom(30, 0.5, 0.1/30)
#' quasibinom_FindParameters(30, 16.6, 9.07)

quasibinom_FindParameters <- function(size, mean, var){
	par <- stats::optim(c(0.5, 0),
		 fn = function(par){
		 	phi <- par[2]*(par[1]/size)
		 	(mean - Equasibinom(size,par[1],phi))^2+
		 		(var - Vquasibinom(size,par[1], phi))^2
		      },
		 lower = c(0,0),
		 upper = c(1,0.5),
		 method = "L-BFGS-B")$par
	par[2] <- par[2]*(par[1]/size)
	names(par) <- c("prob","phi")
	return(par)
}
