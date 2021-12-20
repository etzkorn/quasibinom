#' Quasi-Binomial Expectation and Variance
#'
#' @description Convenience functions for computing the expectation and variance of the quasi-binomial distribution.
#'
#' `Equasibinom` and `Vquasibinom` give the expectation and variance for a given parameter combination
#' `prob`, `phi`, `size`.
#'
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
#' @examples
#'
#' Equasibinom(30, 0.5, 0.1/30)
#' Vquasibinom(30, 0.5, 0.1/30)
#' quasibinom_FindParameters(30, 16.6, 9.07)

Equasibinom <- function(size, prob, phi){
	sum((0:size)*dquasibinom(0:size,size, prob, phi))
}
Vquasibinom <- function(size, prob, phi){
	sum((0:size)^2*dquasibinom(0:size,size, prob, phi)) - Equasibinom(size, prob, phi)^2
}
