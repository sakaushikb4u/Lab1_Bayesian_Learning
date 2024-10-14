
# ===================================================================
# Init Library
# ===================================================================
#knitr::opts_chunk$set(echo = TRUE)
library(ggplot2)
library(gridExtra)

# ===================================================================
# Question 1. Daniel Bernoul 
# ===================================================================

#--------------------------------------------------------------------
# Init code for question 1
#--------------------------------------------------------------------
rm(list = ls())

#--------------------------------------------------------------------
# Code for question 1a
#--------------------------------------------------------------------
# given values
nDraws_1a <- 10000
alpha0 <- 8
beta0 <- 8
s <- 22
n <- 70
f <- n - s
alpha_new <- alpha0 + s
beta_new <- beta0 + f

# calculate the true mean and standard deviation
true_mean <- alpha_new/(alpha_new + beta_new)
true_var <- (alpha_new*beta_new) / ((alpha_new + beta_new)^2 *
            (alpha_new + beta_new +1))
true_sd <- sqrt(true_var)

set.seed(12345)

# init the vectors to store the means and sds
means <- c()
sds <- c()

# calculate means and sds
for (i in 1:nDraws_1a) {
  beta_draws <- rbeta(n = i, shape1 = alpha_new, shape2 = beta_new)
  means[i] <- mean(beta_draws)
  sds[i] <- sd(beta_draws)
}

# plot
df_1a <- data.frame(x = 1:nDraws_1a,
                   true_mean = true_mean,
                   true_sd = true_sd)
# create plot of mean and sd to show the converges status
mean_plot <- ggplot(data = df_1a) + 
  geom_point(aes(x = x, y = means, colour = "Samples"),size=0.1) +
  geom_hline(aes(yintercept = true_mean, colour = "True mean")) + 
  ggtitle("mean converges to the true values") + xlab("nDraws") + ylab("y")
sd_plot <- ggplot(data = df_1a) + 
  geom_point(aes(x = x, y = sds, colour = "Samples"),size=0.1)+
  geom_hline(aes(yintercept = true_sd,colour = "True sd")) +
  ggtitle("sd converges to the true values") + xlab("nDraws") + ylab("y")
grid.arrange(mean_plot, sd_plot, nrow = 2)

#--------------------------------------------------------------------
# Code for question 1b 
#--------------------------------------------------------------------
set.seed(12345)

nDraws_1b <- 10000
sample_b <- rbeta(n = nDraws_1b, shape1 = alpha_new, shape2 = beta_new)
prob_sample_b <- length(sample_b[which(sample_b > 0.3)])/ nDraws_1b
exact_value_1b  <- 1 - pbeta(q = 0.3, shape1 = alpha_new, shape2 = beta_new)

#--------------------------------------------------------------------
# Code for question 1c
#--------------------------------------------------------------------
set.seed(12345)
nDraws_1c <- 10000
sample_1c <- rbeta(n = nDraws_1c, shape1 = alpha_new, shape2 = beta_new)
odds <- sample_1c/(1-sample_1c)

# plot
hist(odds, probability= TRUE, breaks = 100)
lines(density(odds), col = "red")


# ===================================================================
# Question 2. Log-normal distribution and the Gini coefficient
# ===================================================================

#--------------------------------------------------------------------
# Init code for question 2
#--------------------------------------------------------------------
rm(list = ls())

#--------------------------------------------------------------------
# Code for question 2a: calc tau value
#--------------------------------------------------------------------
set.seed(12345)

# given values
observations <- c(33, 24, 48, 32, 55, 74, 23, 17)
n_observations <- length(observations)
nDraws_2a <- 10000
mu <- 3.6
nDraws_2b <- 10000

# calculate tau square value and sigma_square
tau_square <- sum((log(observations)-mu)^2) / n_observations
X <- rchisq(nDraws_2b, n_observations)
sigma_square <- n_observations * tau_square / X


# plot
hist(sigma_square, probability= TRUE, breaks = 100)
lines(density(sigma_square), col = "red")

#--------------------------------------------------------------------
# Code for question 2b: compute the posterior distribution of the Gini
# coefficient
#--------------------------------------------------------------------
phi <- pnorm(sqrt(sigma_square/2), mean = 0, sd = 1)
G <- 2 * phi - 1

# Plot
hist(G, probability= TRUE, breaks = 100)
lines(density(G), col = "red")

#--------------------------------------------------------------------
# Code for question 2c: compute a 95% equal tail credible interval for G
#--------------------------------------------------------------------
credible_interval <- quantile(G, c(0.025, 0.975))

# plot
hist(G, probability = TRUE, breaks = 100)
lines(density(G), col = "red")
abline(v = credible_interval, col = "blue")

#--------------------------------------------------------------------
# Code for question 2d: compute a 95% Highest Posterior
# DensityInterval (HPDI) for G
#--------------------------------------------------------------------
density_G <- density(G)

# find the HPDI
# create the density data frame
density_G_df <- data.frame("x" = density_G$x, "y" = density_G$y)
# order it by y using descending order
density_G_df <- density_G_df[order(density_G_df$y, decreasing = TRUE),]
# calculate the cumulative sum
density_G_df$cumsum <- cumsum(density_G_df$y)
# find the 95% HPDI
HPDI <- density_G_df$x[density_G_df$cumsum < 0.95 * sum(density_G_df$y)]
hdpi_low <- min(HPDI)
hdpi_high <- max(HPDI)

# plot
hist(G, probability= TRUE, breaks = 100)
lines(density(G), col = "red")
abline(v = credible_interval, col = "blue")
abline(v = hdpi_low, col = "green", lty = 2)
abline(v = hdpi_high, col = "green", lty = 2)

# ===================================================================
# 3. Bayesian inference for the concentration parameter in the von
#    Mises distribution
# ===================================================================

#--------------------------------------------------------------------
# Init code for question 3
#--------------------------------------------------------------------
rm(list = ls())

#--------------------------------------------------------------------
# Code for question 3a: calc the posterior
#--------------------------------------------------------------------
# given values
degrees <- c(20, 314, 285, 40, 308, 314, 299, 296, 303, 326)
n_degrees <- length(degrees)
radians <- c(-2.79, 2.33, 1.83, -2.44, 2.23, 2.33, 2.07, 2.02, 2.14, 2.54)
mu_3a <- 2.4
mu_sample <- mean(radians)
sd_sample <- sd(radians)

radians <- (radians-mu_sample)/sd_sample
mu_3a <- (mu_3a-mu_sample)/sd_sample
n_radians <- length(radians)

# grid of k values
k <- seq(from = 0.1,to = 10, by = 0.01) # k>0

# calculate the posterior
result_3a <- exp(k * sum(cos(radians - mu_3a)) - k / 2) /
  besselI(x = k, nu = 0) ^ n_radians
k_max_index <- which.max(result_3a)

#-----------------------------------------
# Code for question 3a: plot the posterior
#-----------------------------------------
# create plot data
df_3a <- data.frame("k" = k, "posterior" = result_3a)

plot_3a <- ggplot(data = df_3a, aes(x=k, y=result_3a)) +
  geom_line() +
  ylab("Posterior value") +
  xlab("k value") +
  ggtitle("Posterior distribution of k")

plot_3a + geom_vline(aes(xintercept = k[k_max_index]))


## (b) Find the (approximate) posterior mode of k from the information in a.

#-----------------------------------------
# Code for question 3b: posterior mode of k 
#-----------------------------------------
k_mode <- k[k_max_index]
