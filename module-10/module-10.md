Module 10
================

Statistical Inference and Basic Hypothesis Testing
==================================================

Preliminaries
-------------

-   Install this package in ***R***: {curl}

Summary
-------

Z and T Statistics for the One Sample Case (testing an observed sample mean for a continuous normally distributed variable compared to a null expectation)

T =

where: {x} = mean of sample observations \_0 = expected mean s = sample standard deviation n = number of sample observations

Z and T Statistics for the Two Sample Case (comparing two means for a continuous normally distributed variable)

Assumptions: - Dealing with normally distributed continuous variables (or those that can be approximated closely by the normal distribution) - When sample size &gt; 30 we can use the Z distribution, but for &lt; 30, use the T distribution

-   CI = mean ± T 1-alpha/2 x SE
-   REJECT Ho if 1-alpha CI around test statistic does not include zero
-   REJECT Ho if p value for obtaining the given test statistic is &lt; alpha

For nonnormally distributed variables...

The same principles apply. Let's consider a proportion...

Or a Poisson CI...

Recall the the Poisson distribution models counts and is determined by a single parameter, *λ*. The Poisson distribution is also useful for modeling rates, or counts that occur over units of time. If we imagine a variable X ∼ Poisson(λt) where λ = E\[X/t\] = the expected count per unit of time and t is the total monitoring time and

Remember that if X ∼ Poisson(λt) then our estimate of λ is λˆ = X/t. Furthermore, we know that V ar(λˆ) = λ/t and so the natural estimate is λˆ/t. While it’s not immediate how the CLT applies in this case, the interval is of the familiar form So our Poisson interval is: Example Estimate ± Z1−α/2SE. λˆ ± Z 1 − α / 2 λˆ

What about other data types, like proportions?

Remember for proportion data, the average proportion = *π* and the standard error of the average proportion was sqrt(*π*(1 − *π*)/*n*)

The Z statistic is equivalent to the T statistic

Errors and Power
----------------

Let's return to the concepts of error and power. Recall that Type I error occurs when you incorrectly reject a true *H*<sub>0</sub>. In any given hypothesis test, the probability of a Type I error is equivalent to the significance level, *α*. Type II error occurs when you incorrectly fail to reject a false *H*<sub>0</sub> (in other words, fail to find evidence in support of a true *H*<sub>*A*</sub>). Since we do not know what the true *H*<sub>*A*</sub> actually is, the probability of committing such an error, labeled *β*, is not usually known in practice.

### Type I Error and the Multiple Testing Problem

Because of how we define *α*, the chance probabilty of falsely rejecting *H*<sub>0</sub>, we would expect to find some "significant" results if we run enough independent hypothesis tests. For example, if we set *α* at 0.05, we expect to find one "significant" result in roughly every 20 tests we run, just by chance. The relation of *α* to the distribution of a variable under a null hypothesis (*μ* = *μ*<sub>0</sub>) versus an alternative hypothesis (*μ* &gt; *μ*<sub>0</sub>) is shown in the figure below (this is an example for a one-tailed test). It should be clear that we can reduce the chance of Type I error by decreasing *α*.

<img src="img/typeI.png" width="600px"/>

Let's explore this via simulation.

We will write some code to simulate a bunch of random datasets from a normal distribution where we know the population mean and standard deviation and then calculate a T statistic and p value for each one. We will then look at the "Type I" error rate... the proportion of times that, based on our sample, we would conclude that it was not drawn from the distribution we know to be true.

First, let's set up a skeleton function we will call `typeI()`. It should take, as arguments, the parameters of the normal distribution we want to simulate from, our sample size, our alpha level, and the number of simulated datasets we want to generate. Type in the code below (and note that we set default values for *α* and number of simulations).

``` r
> typeI <- function(mu0, sigma, n, alpha = 0.05, sims = 10000) {
+ }
```

Now, we will add the body of the function

``` r
> typeI <- function(mu0, sigma, n, alpha = 0.05, k = 10000) {
+     p <- rep(NA, k)  # sets up a vector of empty p values
+     for (i in 1:k) {
+         # sets up a loop to run k simulations
+         x <- rnorm(n = n, mean = mu0, sd = sigma)  # draws a sample from our distribution
+         m <- mean(x)  # calculates the mean
+         s <- sd(x)  # calculates the standard deviation
+         t <- (m - mu0)/(s/sqrt(n))  # calculates the T statistic
+         p[i] <- pt((m - mu0)/(s/sqrt(n)), df = n - 1, lower.tail = FALSE)  # calculates the associated p value
+         # print(p[i])
+     }
+     return(length(p[p < alpha])/k)  # returns the proportion of simulations where p < alpha
+ }
```

Now, run our typeI error test function with a couple of different values of mu0 and sigma and alpha. What error rate is returned? It should be close to alpha!

``` r
> e <- typeI(mu0 = 3, sigma = 4, alpha = 0.05, n = 10)
> e
```

    ## [1] 0.0467

We can address the multiple testing problem by using what is called the Bonferroni correction, which suggests that when doing a total of *N* independent hypothesis tests, each with a significance level of *α*, we should adjusted the *α* level we use to interpret statistical significance as follow: *α*<sub>*B*</sub> = *α*/*N*.

Type II Error
-------------

By reducing the *α* level we use as our criterion for statistical significance, we can reduce the chance of committing Type I, but doing so directly increases our chance of committing a Type II error. The shaded area in this figure below, *β*, is the probability of incorrectly failing to reject the null...

<img src="img/typeII.png" width="600px"/>

It should be clear in this figure that if the critical value (which is defined by *α*) is shifted right or *μ* under the alternative hypothesis shifts left, then *β*, the area under the curve to the left of the critical value, increases! Intuitively, this makes sense: the lower the difference between the true *μ*<sub>*A*</sub> value and *μ*<sub>0</sub> and/or the smaller the *α* level, the harder it will be to reject the null hypothesis that *μ* = *μ*<sub>0</sub>.

Beta usually can’t be calculated in practice because of the need to know what the true distribution actually is. But we can calculate the critical value for our given data (the dashed vertical line) and then use that to calculate what beta would be under a given alternate mean, assuming the same standard error

How does increasing the variability in a given sample impact the Type II error rate? INCRASES it. Larger sample size =&gt; decreased beta

the Type I error rate matches the predefined significance level and so can be decreased by reducing ↵. In contrast, controlling the Type II error rate is a complex balancing act that can involve sample size, significance level, observation variability, and magni- tude of the difference between the true value and the null.

Power

Power is the probability of correctly rejecting a null hypothesis that is untrue. For a test that has a Type II error rate of beta , the statistical power is found simply with 1 - beta. Power of .8 or greater is considered high.

Simulating power...

Let's write a function to do a power simulation for samples of different size. We want a function where we can give it our sample mean, our sample standard deviation, our sample size, and an alternative mean and have it return power.

n &lt;- 1:100

p.test &lt;- function(x){sample.size=n, mu0=mu0, muA=muA,sigma=sigma, alpha=alpha}

typeIIerror &lt;- function(mu0,muA,sigma,n,alpha=0.05,sims=10000){ p &lt;- rep(NA,sims) \# sets up a vector of empty p values for(i in 1:sims){ x &lt;- rnorm(n=n,mean=muA,sd=sigma) m &lt;- mean(x) s &lt;- sd(x) p\[i\] &lt;- 1-pt((m-mu0)/(s/sqrt(n)),df=n-1) \# print(p\[i\]) } return(length(p\[p&gt;=alpha\])/sims) }

mu0 &lt;- 0 muA &lt;- 6 sigma &lt;- 4 n &lt;- 5

typeIIerror(mu0,muA,sigma,n,alpha,sims)

typeIerror(mu0=0,sigma=1,n=40,alpha=0.05) \# simulates 1000 samples of size 40 from random normal distribution
