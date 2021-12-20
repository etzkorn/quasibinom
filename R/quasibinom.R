#' Quasi-Binomial Distribution
#'
#' @description Density, distribution function, quantile function for quasi-binomial distribution.
#'
#' `dquasibinom`, `pquasibinom`, and `qquasibinom` give the density, cumulative density,
#' and quantile distribution functions.
#'
#' @param k vector of quantiles.
#' @param size number of trials (greater than 0).
#' @param prob probability of success on each trial only when phi = 0.
#' @param phi dispersion parameter. Distribution is equivalent to a binomial when phi is 0.
#' Phi must be <= min(prob/size, 1 - prob/size) for a proper density function.
#'
#' @return A numeric vector.
#'
#' @export
#'
#' @details
#' p(x) =  	\choose{\text{size},k} \cdot \text{prob} \cdot
#'          (\text{prob} + x\cdot\phi)^(x-1) \cdot (1-\text{prob}-x\cdot\phi)^(\text{size}-x)
#'
#' The function `quasibinom_FindParameters` uses optim to find the parameter
#' combination to match a given expectation and variance. Answer should not be
#' interpreted as exact, so proceed with caution.
#'
#' @examples
#'
#' dquasibinom(1:30, 30, 0.5, 0.5/30)
#' pquasibinom(1:30, 30, 0.5, 0.5/30)
#' qquasibinom(0.5, 30, 0.5, 0.5/30)

dquasibinom <- function(k, size, prob, phi){
	if(phi > min(prob/size, 1 - prob/size)){
		warning("Phi must be <= min(prob/size, 1 - prob/size) for a proper density function.")
	}
	choose(size,k) * prob * (prob + k*phi)^(k-1) * (1-prob-k*phi)^(size-k)
}
pquasibinom <- function(k, size, prob, phi){
	sapply(k, function(x) sum(dquasibinom(0:x,size, prob, phi)))
}
qquasibinom <- function(p, size, prob, phi){
	min((0:size)[pquasibinom(0:size, size, prob, phi) >= p ])
}
