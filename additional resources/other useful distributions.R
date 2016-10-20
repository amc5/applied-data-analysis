module-10

Other Useful Distributions


Negative Binomial (discrete)

Similar to the binomial and poisson, this distribution is derived from the number of failures rather than successes. The distribution maps the number of Bernoulli trials until a specified event occurs r times. The negative binomial differs from the poisson in that the variance may be different from the expected value.

Examples include the number of heads before one gets 2 tails; number of individuals per patch; number of successful trapping events in a trap day. In ecology the negative binomial is useful for modeling processes where the variance is expected to be larger than the mean (overdispersed) and there need not be an upper limit. The negative binomial, like the binomial, is bounded at 0 (no negative values). The negative binomial arrises in situations similar to the poisson but where lambda is itself variable.

The negative binomial is represented in one of two ways

X ~ Nbinom(p,n)

Where p is the probability per trial and n is the number of successes awaited. Or:

	X ~ Nbinom(μ, k)

Where μ is the mean number of failures expected and k is the overdispersion parameter, which measures the amount of clustering or aggregation or heterogeneity.

For additional details see (Bolker, 2008)

Note that there are two conflicting meanings for overdispersion. In the statistical sense, a distribution such as the negative binomial may be called overdispersed when its variance is greater than the mean. As discussed above for the poisson distribution, statistical overdispersion occurs when events are clustered or aggregated.

In spatial statistics, overdispersed means that occurrences are less aggregated, more evenly distributed on the landscape than one would expect by chance.

The R functions for this distribution are: [p,d,q,r]nbinom(), which are included in the base package.

Geometric (discrete)

The number of Bernoulli trials until a specified event occurs. The goemetric distribution is a special case of the negative binomial, with k or n=1.

Ecological examples include number of primates sampled until we get one that is parasitized; waiting time until a mutation occurs at a nucleotide; waiting time until a female accepts a mate; survival analysis of prey when a predator is present.

R has the following functions for this distribution as part of the base package: [p,d,q,r]geom()

Hypergeometric

Like the binomial, but sampling from a small population without replacement, so each success changes the probability of success in the next trial. This is the basis for mark-recapture studies of population size.

R has [d,p,q,r]hyper()

There is also a generalized hypergeometric function [d,p,q,r]ghyper() in the SuppDists package.

Multinomial (discrete)

Like a binomial, but we are counting more than two possible outcomes. With c different possible outcomes or categories {c1, c2, …, cc}, we are counting the number of each result n1, n2, …., nc . The probability of each outcome is {π1, π2, … , πc} and formula. For n independent observations (trials) where formula the probability of there being n1 in c1 is formula

R has two functions for the multinomial distribution rmultinom dmulitnom

The multnomial is a multivariate distribution and hence there is no pmultinom or qmultinom because there are multiple ways to increase values in any category.

The multinomial can be used to model species relative abundance.

Example: Simulate a sample of 100 individuals taken from a community of 5 species with the following relative abundances, {0.6, 0.2, 0.1, 0.06, 0.04}

tprobs <- c(0.6, 0.2, 0.1, 0.06, 0.04)
sum(tprobs)  # must equal 1
## [1] 1
tdist <- rmultinom(1, 100, tprobs)
print(tdist)
##      [,1]
## [1,]   59
## [2,]   20
## [3,]   15
## [4,]    5
## [5,]    1





Gamma (continuous)

The Gamma distribution is the distribution of waiting times until a certain number of events occur. It is the continuous version of the negative binomial. Examples include the number of days till a certain number of individuals in a population die. Gamma distributions are described by two variables: α = scale (the length of an event) and β = shape (number of events). The scale can alternately be given as a rate (1/scale). R base provides [d,p,q,r]gamma().

Exponential Distribution (continuous)

A special case of the Gamma distribution where shape = 1. Often used to describe the waiting time until an event occurs (so, this is the continuous analog to the Geometric). Used to describe radioactive decay, extinction times, durations of species in the fossil record, persistence of neutral mutations in a population.

A key feature is that the exponential distribution is “memoryless”, that is the waiting time to the next event is totally independent of how long it has been since the last one. This is a classic statistical problem in gambling, where many people incorrectly think that a number which hasn’t been hit in roulette for a long time should come up sooner.

R provides [d,p,q,r]exp()

Example 1: If species longevity averages 2 Ma in primates (i.e. extinction rate is 0.5 per Ma), generate a random sample of 20 species duration periods in Ma.

rexp(100, rate = 0.5)
##   [1] 0.1670 1.2010 0.2761 1.2009 5.3202 1.4862 0.2420 1.1417 3.4528 2.3997
##  [11] 0.3255 0.6811 0.1693 1.4625 1.3807 1.9729 2.8597 3.8163 0.9027 0.0666
##  [21] 1.2326 0.5042 0.2857 1.8312 1.4769 0.1040 0.1665 1.7412 0.6210 1.1869
##  [31] 4.0623 2.0510 0.4389 0.9389 0.7407 1.2196 3.3009 1.0369 2.1738 4.1663
##  [41] 1.7728 4.0732 0.6512 0.5474 1.2495 2.6952 0.1115 0.6233 1.1816 1.5437
##  [51] 6.0392 0.8166 5.9265 4.2570 2.6502 3.1559 5.4457 1.3616 8.7826 1.5311
##  [61] 2.2538 4.9279 0.4289 2.5200 1.0100 1.3147 1.0356 0.3052 3.3199 1.2784
##  [71] 2.0871 0.5249 2.5888 0.4483 1.5446 0.8165 0.7201 0.8123 0.3498 2.0232
##  [81] 5.8709 2.8607 2.2707 1.3702 1.2338 0.0109 1.9479 1.2180 3.8666 0.6308
##  [91] 0.7173 0.8966 0.6080 4.4530 7.9656 0.1614 0.2933 0.1730 0.7527 2.3539
Example 2: Plot the probabilities of various species durations in Ma if the average extinction rate is 0.5 spp. Ma-1. The second figure shows the cumulative probabilities. 95% of all species will be extinct after ca 6 Ma. .

plot(seq(0, 100, 0.1), dexp(seq(0, 100, 0.1), rate = 0.5), type = "l")
plot of chunk exp_density

qexp(0.95, 0.5)
## [1] 5.991
plot(seq(0, 100, 0.1), pexp(seq(0, 100, 0.1), rate = 0.1), type = "l")
plot of chunk exponential2

Chi Square

The chi-square distribution has a single parameter ν called the degrees of freedom. The chi-square distribution is a special case of the Gamma where the shape is 2 and the scale is half the degrees of freedom (i.e. ν/2). This distribution arises often in statistical inference and is closely tied to the normal distribution. R has functions [d,p,q,r]chisq()

Log Normal (continuous)

Taking the log of a lognormal r.v. produces a normal r.v. The log normal arises as the product of many independent random variables, similar to the way the normal distribution results from the summation of many independent random variables.



Here are a couple of important things to remember...

### Discrete Random Variables

| Name           | Notation                | PMF $f(x)$ = $P(X=x)$                      | $\mu$ = $E(X)$     | $\sigma^2$ = *Var* $(X)$       |                          |
	|:-------------- |:----------------------- |:------------------------------------------ |:------------------ |:------------------------------ |:-----------------------  |
	| Bernoulli      | $X$ ~ *BERN*$(p)$       | $f(x)$ = $p^x(1-p)^{1-x}$                  | *p*                | *p(1-p)*                       | x = {0,1}                |
	| Binomial       | $X$ ~ *BIN*$(n, p)$     | <img src="img/binom-1.svg" width="220px"/> | *np*               | *np(1-p)*                      | x = {0,1,..., k}         |
	| Poisson        | $X$ ~ *POIS*$(\lambda)$ | <img src="img/poisson.svg" width="175px"/> | $\lambda$          | $\lambda$                      | x = {0,1,..., +$\infty$} |

	### Continuous Random Variables

	| Name         | Notation                      | PDF $f_X(x)$                                 | $\mu$ = $E(X)$     | $\sigma^2$ = *Var* $(X)$ |
	|:------------ |:----------------------------- | -------------------------------------------- |:------------------ |:------------------- |:---------------------- |
	| Beta         | $X$ ~ *BETA*$(\alpha,\beta)$  | $f(x)$ = $K$ $x^{\alpha-1}(1-x)^{\beta-1}$   |                    |                     |                        |
	| Uniform      | $X$ ~ *U*$(a,b)$              | <img src="img/uni-1.svg" width="110px"/>     | <img src="img/uni-2.svg" width="100px"/> | <img src="img/uni-3.svg" width="115px"/> |                        |
	| Normal      | $X$ ~ *N*$(\mu,\sigma^2)$      |                                              |                    |                     |                        |
	| Exponential |
	| Student's t |


Testing for normality statistically...


The Shapiro-Wilk Test yields a test statistic W, by finding the largest deviation from the expected line in a qqplot. This is a powerful test for normality, but does not work well when there are many ties in the data. This test is not part of the base installation of R and requires the cwhmisc package.

library(cwhmisc)

shapiro.test(x)

Another powerful test for normality is the Anderson-Darling test. This test comparable in power to the SW test. It requires the nortest package.

library(nortest)

Testing for fit to a distribution of your choice. A more general test is the Kolmogorov-Smirnov Test (KS test). For details, see Zar, Biostatistical Analysis

KS Test. If you test for a normal distribution, it assumes you are dealing with a standard normal distribution (mean 0, standard deviation 1).

You can also use ks.test to compare 2 distributions to determine whether they come from the same distribution: warning though, if means or variances are different, the ks test will say they are different distributions.

ks.test(x, rnorm(100, 10, 3))

The K-S Test can be applied for most distributions. Its main limitation is that real statisticians will tell you that you should not apply the K-S test using a distribution with parameters estimated directly from the data. That is, you shouldn’t estimate the mean and variance of a binomial distribution, and use that to generate the null distribution with which you compare your data. In practice, many biologists ignore this stricture, because the KS test is more powerful than the alternative, the chi squared test.

The KS test can be used for almost any distribution in R, you just need to replace “pnorm” with the distribution you want.



Coefficient of Variation = ratio of the standard deviation to the mean

v <- rnorm(100, mean=5, sd=2)
CV <- sd(v)/mean(v)


#### The Chi Squared Distribution

The **Chi Squared Distribution** models sums of squared normal variates and is thus often used in tests concerning sample variances of normally distributed data. Like the t distribution, it is dependent upon specification of a degrees of freedom using the `df=` argument.
