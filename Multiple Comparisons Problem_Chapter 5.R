#### This R-script is for the paper 'Multiple Comparisons Problem:  
#### Bonferroni Correction vs. Benjamini-Hochberg Procedure', Chapter 5.
### variable legend:
# m ... number of hypotheses testings
# p ... vector of p-values
# alpha ... significance level for hypothesis testing
# X ... matrix of simulated data
# n ... number of variables
# r ... number of runs/repetitions
# H0 ... number of true null hypotheses
# H1 ... number of true alternative hypotheses


####################### (5.1) Application of both Procedures: ##################
rm(list=ls())
set.seed(36)


r <- 2500
n <- 36
m <- n # in this setting n = m
alpha <- 0.05
p <- matrix(NA, nrow = r, ncol = n)

# test simultaneously if mean of all variables equal 0
# H0_i: mu_1 = ... = mu_n = 0 with i in 1:m
# if one H0_ij is rejected, so is H0_i (j in 1:n)
for(i in 1:r){
  # normally distributed random deviates
  X <- replicate(n, expr = rnorm(1000, 0, 1))  # design matrix X
  # test per column, resulting in m p-values per run (r runs total)
  p[i, ] <- apply(X = X, MARGIN = 2, 
                  FUN = function(x){
                    t.test(x = x, mu = 0, conf.level = 1-alpha)$p.value
                    })
}

# in m runs, how often did we reject (expected: 5% of time = Type-1 error)
bool <- rep(NA, r)

for(i in 1:r){
  bool[i] <- any(p[i, ] < 0.05) 
}                      
mean(bool)
sum(bool) # = r*mean(bool)


## Bonferroni correction
bonf_bool <- rep(NA, r)

for(i in 1:r){
  bonf_bool[i] <- any(p[i, ] < 0.05/m) 
}                      
mean(bonf_bool)
sum(bonf_bool)


## Benjamini-Hochberg procedure
p_ord <- matrix(NA, nrow = r, ncol = n)
crit_val <- matrix(NA, nrow = r, ncol = n)

for(i in 1:r){
  # order p-values in each row
  p_ord[i, ] <- p[i, ][order(p[i, ])]
  # calculate critical value, i.e. (j*alpha)/m
  for(j in 1:n){
    crit_val[i, j] <- (j * alpha)/m  
  }
}

# how many rejected
bh_bool <- rep(NA, r)

for(i in 1:r){
  bh_bool[i] <- any(p_ord[i, ] < crit_val[i, ]) 
}                      
mean(bh_bool)
sum(bh_bool)






