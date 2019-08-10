library(pcev)
library(RMTstat)

set.seed(12345)
p = 20; n = 100

Y <- matrix(rnorm(n*p), nrow = n)
X <- rnorm(n)

pcev_out <- computePCEV(Y, X)

p = 200
Y <- matrix(rnorm(n*p), nrow = n)

pcev_out2 <- computePCEV(Y, X, estimation = "singular")

# Estimating the null distribution----
null_lambdas <- replicate(50, {
  X_perm <- sample(X)
  computePCEV(Y, X_perm, 
              estimation = "singular", 
              nperm = 1)$largestRoot
})

hist(log(null_lambdas), freq = FALSE)

# Method of moments
sample_mean <- mean(log(null_lambdas))
sample_sd <- sd(log(null_lambdas))

theo_mean <- -1.2065335745820
theo_sd <- sqrt(1.607781034581)

sigma <- sample_sd/theo_sd
mu <- sample_mean - sigma * theo_mean

transformed_lambdas <- (log(null_lambdas) - mu)/sigma

hist(transformed_lambdas,
     freq = FALSE)

x <- seq(-5, 2, length=1000)
y <- dtw(x)
lines(x, y, col = 'blue')

plot(ecdf(transformed_lambdas),
     main = "")
lines(x, ptw(x), col = 'blue')

# BLK gene----
library(ggplot2)
library(tidyr)

BLK_boundaries <- c(11235000, 11385000)/1000000

output <- computePCEV(t(methylation2), pheno2,
                      estimation = "singular",
                      nperm = 100)

data_plot <- data.frame(pos = position2/1000000,
                        PCEV = output$VIMP)

ggplot(data_plot, aes(pos, PCEV)) + 
  geom_point() + 
  theme_minimal() +
  ylim(c(0, 1))

# Compare to PCA
pca_out <- prcomp(t(methylation2))

data_plot$PCA <- abs(pca_out$rotation[,1])
data_plot$PCA <- data_plot$PCA *(max(data_plot$PCEV)/max(data_plot$PCA))

gather(data_plot, key = "Model",
       value = "VIP", PCEV, PCA) %>% 
  ggplot(aes(pos, VIP)) + 
  geom_point() + 
  geom_segment(y = 0.9, yend = 0.9,
               x = BLK_boundaries[1],
               xend = BLK_boundaries[2],
               colour = 'blue') +
  theme_minimal() +
  ylim(c(0, 1)) +
  facet_grid(~Model)

