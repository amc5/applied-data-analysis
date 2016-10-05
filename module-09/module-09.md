Module 09
================

Statistical Inference and Basic Hypothesis Testing
==================================================

Preliminaries
-------------

-   Install this packages in ***R***: {manipulate}

Objectives
----------

> The objective of this module is to begin our discussion of *statistical inference* from a frequentist/classical statistics approach. Doing so means that we need to cover basics of probability and distributions.

Review of Common Distributions
------------------------------

Here are a couple of important things to remember...

### Discrete Random Variables

| Name      | Notation              | PMF *f*(*x*) = *P*(*X* = *x*)                           | *μ* = *E*(*X*) | *σ*<sup>2</sup> = *Var* (*X*) |                   |
|:----------|:----------------------|:--------------------------------------------------------|:---------------|:------------------------------|:------------------|
| Bernoulli | *X* ~ *BERN*(*p*)     | *f*(*x*) = *p*<sup>*x*</sup>(1 − *p*)<sup>1 − *x*</sup> | *p*            | *p(1-p)*                      | x = {0,1}         |
| Binomial  | *X* ~ *BIN*(*n*, *p*) | <img src="img/binom-1.svg" width="220px"/>              | *np*           | *np(1-p)*                     | x = {0,1,..., k}  |
| Poisson   | *X* ~ *POIS*(*λ*)     | <img src="img/poisson.svg" width="175px"/>              | *λ*            | *λ*                           | x = {0,1,..., +∞} |

### Continuous Random Variables

| Name        | Notation                        | PDF *f*<sub>*X*</sub>(*x*)                                      | *μ* = *E*(*X*)                           | *σ*<sup>2</sup> = *Var* (*X*)            |     |
|:------------|:--------------------------------|-----------------------------------------------------------------|:-----------------------------------------|:-----------------------------------------|:----|
| Beta        | *X* ~ *BETA*(*α*, *β*)          | *f*(*x*) = *K* *x*<sup>*α* − 1</sup>(1 − *x*)<sup>*β* − 1</sup> |                                          |                                          |     |
| Uniform     | *X* ~ *U*(*a*, *b*)             | <img src="img/uni-1.svg" width="110px"/>                        | <img src="img/uni-2.svg" width="100px"/> | <img src="img/uni-3.svg" width="115px"/> |     |
| Normal      | *X* ~ *N*(*μ*, *σ*<sup>2</sup>) |                                                                 |                                          |                                          |     |
| Exponential |                                 |
| Student's t |                                 |

Coefficient of Variation = ratio of the standard deviation to the mean

v &lt;- rnorm(100, mean=5, sd=2) CV &lt;- sd(v)/mean(v)

Roughly, the central limit theorem states that the distribution of the sum (or average) of a large number of independent, identically distributed variables will be approximately normal, regardless of the underlying distribution.

The central limit theorem states that averages of random variables independently drawn from independent distributions converge in distribution to the normal, that is, become normally distributed when the number of random variables is sufficiently large. Also, variables that are expected to be the sum of many independent processes (such as measurement errors) often have distributions that are nearly normal.
