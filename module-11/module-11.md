Module 11
================

Type I and II Errors and Statistical Power
==========================================

Preliminaries
-------------

-   Install these packages in ***R***: {curl}, {ggplot2}, {manipulate}

Overview
--------

Let's return to the concepts of **error** and **power**. Recall that Type I error occurs when you incorrectly reject a true *H*<sub>0</sub>. In any given hypothesis test, the probability of a Type I error is equivalent to the significance level, *α*, and it is this type of error we are often trying to minimize when we are doing classical statistical inference. Type II error occurs when you incorrectly fail to reject a false *H*<sub>0</sub> (in other words, fail to find evidence in support of a true *H*<sub>*A*</sub>). Since we do not know what the true *H*<sub>*A*</sub> actually is, the probability of committing such an error, labeled *β*, is not usually known in practice.

| What is True      | What We Decide    | Result                                    |
|-------------------|-------------------|-------------------------------------------|
| *H*<sub>0</sub>   | *H*<sub>0</sub>   | Correctly 'accept' the null               |
| *H*<sub>0</sub>   | *H*<sub>*A*</sub> | Falsely reject the null (Type I error)    |
| *H*<sub>*A*</sub> | *H*<sub>0</sub>   | Falsely 'accept' the null (Type II error) |
| *H*<sub>*A*</sub> | *H*<sub>*A*</sub> | Correctly reject the null                 |

Type I Error and the Multiple Testing Problem
---------------------------------------------

Because of how we define *α*, the chance probabilty of falsely rejecting *H*<sub>0</sub> when *H*<sub>0</sub> is actually true, we would expect to find some "significant" results if we run enough independent hypothesis tests. For example, if we set *α* at 0.05, we would expect to find one "significant" result in roughly every 20 tests we run, just by chance. The relation of *α* to the distribution of a variable under a null hypothesis (*μ* = *μ*<sub>0</sub>) versus an alternative hypothesis (e.g., *μ* &gt; *μ*<sub>0</sub>) is shown in the figure below (this is an example for a one-tailed test). It should be clear that we can reduce the chance of Type I error by decreasing *α*.

<img src="img/typeI.png" width="600px"/>

Let's explore this via simulation.

We will write some code to simulate a bunch of random datasets from a normal distribution where we set the expected population mean (*μ*) and standard deviation (*σ*) and then calculate a T statistic and p value for each one. We will then look at the "Type I" error rate... the proportion of times that, based on our sample, we would conclude that it was not drawn from the distribution we know to be true.

First, let's set up a skeleton function we will call `typeI()` to evaluate the Type I error rate. It should take, as arguments, the parameters of the normal distribution for the null hypothesis we want to simulate from (*μ*<sub>0</sub> and *σ*), our sample size, our alpha level, what type of t test we want to do ("greater", "less", "two.tailed") and the number of simulated datasets we want to generate. Type in the code below (and note that we set default values for *α* and the number of simulations).

``` r
> typeI <- function(mu0, sigma, n, type = "two-tailed", alpha = 0.05, k = 10000) {
+ }
```

Now, we will add the body of the function.

``` r
> typeI <- function(mu0, sigma, n, type = "two-tailed", alpha = 0.05, k = 1000) {
+     p <- rep(NA, k)  # sets up a vector of empty p values
+     for (i in 1:k) {
+         # sets up a loop to run k simulations
+         x <- rnorm(n = n, mean = mu0, sd = sigma)  # draws a sample from our distribution
+         m <- mean(x)  # calculates the mean
+         s <- sd(x)  # calculates the standard deviation
+         t <- (m - mu0)/(s/sqrt(n))  # calculates the T statistic for the sample drawn from the null distribution relative to the null distribution
+         if (type == "less") {
+             p[i] <- pt(t, df = n - 1, lower.tail = TRUE)  # calculates the associated p value
+         }
+         if (type == "greater") {
+             p[i] <- pt(t, df = n - 1, lower.tail = FALSE)
+         }
+         if (type == "two-tailed") {
+             if (t > 0) {
+                 p[i] <- 2 * pt(t, df = n - 1, lower.tail = FALSE)
+             }
+             if (t < 0) {
+                 p[i] <- 2 * pt(t, df = n - 1, lower.tail = TRUE)
+             }
+         }
+     }
+     curve(dnorm(x, mu0, sigma), mu0 - 4 * sigma, mu0 + 4 * sigma, main = paste("Type I = ", 
+         length(p[p < alpha])/k, sep = ""), xlab = "x", ylab = "Pr(x)")
+     if (type == "less") {
+         polygon(cbind(c(mu0 - 4 * sigma, seq(from = mu0 - 4 * sigma, to = mu0 - 
+             qnorm(1 - alpha) * sigma, length.out = 100), mu0 - qnorm(1 - alpha) * 
+             sigma)), c(0, dnorm(seq(from = mu0 - 4 * sigma, to = mu0 - qnorm(1 - 
+             alpha) * sigma, length.out = 100), mean = mu0, sd = sigma), 0), 
+             border = "black", col = "gray")
+     }
+     if (type == "greater") {
+         polygon(cbind(c(mu0 + qnorm(1 - alpha) * sigma, seq(from = mu0 + qnorm(1 - 
+             alpha) * sigma, to = mu0 + 4 * sigma, length.out = 100), mu0 + 4 * 
+             sigma)), c(0, dnorm(seq(from = mu0 + qnorm(1 - alpha) * sigma, to = mu0 + 
+             4 * sigma, length.out = 100), mean = mu0, sd = sigma), 0), border = "black", 
+             col = "gray")
+     }
+     if (type == "two-tailed") {
+         polygon(cbind(c(mu0 - 4 * sigma, seq(from = mu0 - 4 * sigma, to = mu0 - 
+             qnorm(1 - alpha/2) * sigma, length.out = 100), mu0 - qnorm(1 - alpha/2) * 
+             sigma)), c(0, dnorm(seq(from = mu0 - 4 * sigma, to = mu0 - qnorm(1 - 
+             alpha/2) * sigma, length.out = 100), mean = mu0, sd = sigma), 0), 
+             border = "black", col = "gray")
+         polygon(cbind(c(mu0 + qnorm(1 - alpha/2) * sigma, seq(from = mu0 + qnorm(1 - 
+             alpha/2) * sigma, to = mu0 + 4 * sigma, length.out = 100), mu0 + 
+             4 * sigma)), c(0, dnorm(seq(from = mu0 + qnorm(1 - alpha/2) * sigma, 
+             to = mu0 + 4 * sigma, length.out = 100), mean = mu0, sd = sigma), 
+             0), border = "black", col = "gray")
+     }
+     return(length(p[p < alpha])/k)  # returns the proportion of simulations where p < alpha
+ }
```

Can you explain what each step of this code is doing?

Now, run our Type I error test function with a couple of different values of *μ*<sub>0</sub>, *σ*, and *α*. What error rates are returned? They should be always be close to *α*!

``` r
> e <- typeI(mu0 = -3, sigma = 2, n = 500, type = "greater", alpha = 0.05)
```

![](img/unnamed-chunk-3-1.png)

``` r
> e
```

    ## [1] 0.044

#### CHALLENGE:

How does the Type I error rate change with *n*?

### Bonferroni Correction

We can address the multiple testing problem by using what is called the Bonferroni correction, which suggests that when doing a total of *N* independent hypothesis tests, each with a significance level of *α*, we should adjusted the *α* level we use to interpret statistical significance as follow: *α*<sub>*B*</sub> = *α*/*N*.

For example, if we run 10 independent hypothesis tests, then we should set our adjusted *α* level for each test as 0.05/10 = 0.005. Note that many statisticians consider the Bonferroni correction to be a particularly conservative one, and there are other corrections we might use to account for multiple testing.

Type II Error
-------------

By reducing the *α* level we use as our criterion for statistical significance, we can reduce the chance of committing Type I, but doing so directly increases our chance of committing a Type II error. The shaded area in this figure below, *β*, is the probability of incorrectly failing to reject the null...

<img src="img/typeII.png" width="600px"/>

It should be clear in this figure that if the critical value (which is defined by *α*) is shifted right or if *μ* under the alternative hypothesis shifts left, then *β*, the area under the curve to the left of the critical value, increases! Intuitively, this makes sense: the lower the difference between the true *μ*<sub>*A*</sub> value and *μ*<sub>0</sub> and/or the smaller the *α* level, the harder it will be to reject the null hypothesis that *μ* = *μ*<sub>0</sub>.

In practice, we cannot usually calculate *β* because of the need to know where the true distribution is really centered (i.e., what is the value of *μ*<sub>*A*</sub>). However, using our data, we can explore what *β* is expected to look like under different sample sizes, *α* levels, and expected *μ*<sub>*A*</sub>.

Let's do this using the simulation approach we developed above. Again, we will write some code to simulate a bunch of random datasets, this time drawn from a normal distribution associated with a particular alternative hypothesis, *H*<sub>*A*</sub>, i.e., where we know that the expected population mean is *μ*<sub>*A*</sub> and the standard deviation is *σ*. We then calculate a T statistic based on each sample dataset relative to *μ*<sub>0</sub>, the expected mean under *H*<sub>0</sub>, and determine the associated p value for each one. Based on this, we can calculate the Type II error rate... the proportion of times that, based on our sample, we would conclude that it was drawn from the *H*<sub>0</sub> distribution rather than the *H*<sub>*A*</sub> distribution that we "know"" to be true.

``` r
> typeII <- function(mu0, muA, sigma, n, type = "two-tailed", alpha = 0.05, k = 1000) {
+     p <- rep(NA, k)  # sets up a vector of empty p values
+     for (i in 1:k) {
+         x <- rnorm(n = n, mean = muA, sd = sigma)  # draw from Ha
+         m <- mean(x)
+         s <- sd(x)
+         t <- (m - mu0)/(s/sqrt(n))  # calculates the T statistic for the sample drawn from Ha relative to the null distribution
+         if (type == "less") {
+             p[i] <- pt(t, df = n - 1, lower.tail = TRUE)  # calculates the associated p value
+         }
+         if (type == "greater") {
+             p[i] <- pt(t, df = n - 1, lower.tail = FALSE)
+         }
+         if (type == "two.tailed") {
+             if (t > 0) {
+                 p[i] <- 2 * pt(t, df = n - 1, lower.tail = FALSE)
+             }
+             if (t < 0) {
+                 p[i] <- 2 * pt(t, df = n - 1, lower.tail = TRUE)
+             }
+         }
+     }
+     
+     curve(dnorm(x, mu0, sigma), mu0 - 4 * sigma, mu0 + 4 * sigma, main = paste("Type II - Prob of Incorrectly Failing To Reject Null\nBeta = ", 
+         length(p[p >= alpha])/k, sep = ""), xlab = "x", ylab = "Pr(x)", col = "blue", 
+         lty = 3, xlim = c(min(c(mu0 - 4 * sigma, muA - 4 * sigma)), max(c(mu0 + 
+             4 * sigma, muA + 4 * sigma))), ylim = c(0, max(c(dnorm(mu0, mu0, 
+             sigma)), dnorm(muA, muA, sigma), dnorm(muA, muA, sigma/sqrt(n)))))
+     curve(dnorm(x, muA, sigma), muA - 4 * sigma, muA + 4 * sigma, col = "red", 
+         lty = 3, add = TRUE)
+     
+     curve(dnorm(x, mu0, sigma/sqrt(n)), mu0 - 4 * sigma/sqrt(n), mu0 + 4 * sigma/sqrt(n), 
+         col = "blue", add = TRUE)
+     curve(dnorm(x, muA, sigma/sqrt(n)), muA - 4 * sigma/sqrt(n), muA + 4 * sigma/sqrt(n), 
+         col = "red", add = TRUE)
+     abline(h = 0)
+     
+     if (type == "less") {
+         polygon(cbind(c(mu0 - qnorm(1 - alpha) * sigma/sqrt(n), seq(from = mu0 - 
+             qnorm(1 - alpha) * sigma/sqrt(n), to = muA + 4 * sigma/sqrt(n), 
+             length.out = 100), muA + 4 * sigma/sqrt(n))), c(0, dnorm(seq(mu0 - 
+             qnorm(1 - alpha) * sigma/sqrt(n), to = muA + 4 * sigma/sqrt(n), 
+             length.out = 100), mean = muA, sd = sigma/sqrt(n)), 0), border = "red", 
+             col = "red")
+         abline(v = mu0 - qnorm(1 - alpha) * sigma/sqrt(n), col = "black", lty = 3, 
+             lwd = 2)
+     }
+     
+     if (type == "greater") {
+         polygon(cbind(c(muA - 4 * sigma/sqrt(n), seq(from = muA - 4 * sigma/sqrt(n), 
+             to = mu0 + qnorm(1 - alpha) * sigma/sqrt(n), length.out = 100), 
+             mu0 + qnorm(1 - alpha) * sigma/sqrt(n))), c(0, dnorm(seq(from = muA - 
+             4 * sigma/sqrt(n), to = mu0 + qnorm(1 - alpha) * sigma/sqrt(n), 
+             length.out = 100), mean = muA, sd = sigma/sqrt(n)), 0), border = "red", 
+             col = "red")
+         abline(v = mu0 + qnorm(1 - alpha) * sigma/sqrt(n), col = "black", lty = 3, 
+             lwd = 2)
+     }
+     
+     return(length(p[p >= alpha])/k)
+ }
> 
> e <- typeII(mu0 = 2, muA = 5, sigma = 3, n = 6, type = "greater")  # Ha > H0
```

![](img/unnamed-chunk-4-1.png)

``` r
> e
```

    ## [1] 0.338

``` r
> e <- typeII(mu0 = 2, muA = 1, sigma = 3, n = 6, type = "less")  # Ha < H0
```

![](img/unnamed-chunk-4-2.png)

``` r
> e
```

    ## [1] 0.826

Power
-----

Power is the probability of correctly rejecting a null hypothesis that is untrue. For a test that has a Type II error rate of *β*, the statistical power is defined, simply, as 1 − *β*. Power values of 0.8 or greater are considered high.

#### CHALLENGE:

Let's explore the effects of increasing the variability (*σ*) in a given sample and the impact of increasing sample size (*n*). How do each of these impact the Type II error rate? How do they effect power?

``` r
> library(ggplot2)
> library(manipulate)
> 
> n <- 1:100
> s <- seq(from = 1, to = 10, by = 1)
> p <- matrix(data = NA, nrow = length(n), ncol = length(s))
> p <- as.data.frame(p)
> for (i in 1:length(n)) {
+     for (j in 1:length(s)) {
+         p[i, j] <- typeII(mu0 = 2, muA = 4, sigma = j, type = "two-tailed", 
+             alpha = 0.05, n = i)
+     }
+ }
> p <- cbind(p, n)
> manipulate(ggplot(data = p) + xlab("sample size n") + ylab("Type II Error Rate (Blue)\nand\nPower (Red)") + 
+     ylim(0, 1) + geom_point(aes(x = n, y = p[, sigma]), colour = "blue", alpha = 1/2) + 
+     geom_line(aes(x = n, y = p[, sigma]), colour = "blue", alpha = 1/2) + geom_line(aes(x = n, 
+     y = 0.8), colour = "red", lty = 3) + geom_point(aes(x = n, y = 1 - p[, sigma]), 
+     colour = "red", alpha = 1/2) + geom_line(aes(x = n, y = 1 - p[, sigma]), 
+     colour = "red", alpha = 1/2), sigma = slider(1, 10, initial = 3, step = 1))
```
