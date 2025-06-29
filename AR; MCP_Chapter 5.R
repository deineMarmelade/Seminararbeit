#### This R-script is for the paper 'Multiple Comparisons Problem:  
#### Bonferroni Correction vs. Benjamini-Hochberg Procedure', Chapter 5.
### variable legend:
# m ... number of hypotheses testings
# p ... matrix of p-values
# alpha ... significance level for hypothesis testing
# X ... matrix of simulated data
# n ... number of variables (number of groups), equals m in this setting 
# r ... number of runs/repetitions


####################### (5.1) Application of both Procedures: ##################
rm(list=ls())
set.seed(36)


r <- 50
m <- n <- 500
alpha <- 0.05
p <- matrix(NA, nrow = r, ncol = n) # rows: repetitions, columns: each group 

# test simultaneously if mean of all groups equal 0
# H0_i: \theta_i = 0 with i in 1:m
# if one H0_i is rejected, so is global H0 per run
for(i in 1:r){
  # normally distributed random deviates
  X <- replicate(n, expr = rnorm(1000, 0, 1))  # design matrix X
  # test per column (group), resulting in m p-values per run (r runs total)
  p[i, ] <- apply(X = X, MARGIN = 2, 
                  FUN = function(x){
                    t.test(x = x, mu = 0, conf.level = 1-alpha)$p.value
                  })
}

# for m groups, how often did we reject (expected: 5% of time = Type-1 error)
lokalH0_bool <- rep(NA, r) # count of rejected H0_i out of 5000 groups per run  
globH0_bool <- rep(NA, r) # count of rejected global H0 per run  

for(i in 1:r){
  lokalH0_bool[i] <- sum((p[i, ] < alpha)) 
  globH0_bool[i] <- any(p[i, ] < alpha)
}                      
mean(lokalH0_bool)
mean(globH0_bool)
sum(lokalH0_bool)
sum(globH0_bool) # = r*mean(globH0_bool)


# what we wanted: 5% of times reject purely by chance
ACTUALlokalH0 <- m*alpha
ACTUALglobH0 <- r*alpha
ACTUALlokalH0
ACTUALglobH0


## Bonferroni correction
bonflokalH0_bool <- rep(NA, r)
bonfglobH0_bool <- rep(NA, r)

for(i in 1:r){
  bonflokalH0_bool[i] <- sum(p[i, ] < alpha/m)
  bonfglobH0_bool[i] <- any(p[i, ] < alpha/m)
}                      
mean(bonflokalH0_bool)
mean(bonfglobH0_bool)
sum(bonflokalH0_bool)
sum(bonfglobH0_bool)


## Benjamini-Hochberg procedure
p_ord <- matrix(NA, nrow = r, ncol = n)
thresh <- matrix(NA, nrow = r, ncol = n)

for(i in 1:r){
  # order p-values in each row
  p_ord[i, ] <- p[i, ][order(p[i, ])]
  # calculate threshold, i.e. (j*alpha)/m
  for(j in 1:n){
    thresh[i, j] <- (j * alpha)/m  
  }
}

# how many rejected
bhlokalH0_bool <- rep(NA, r)
bhglobH0_bool <- rep(NA, r)

for(i in 1:r){
  bhlokalH0_bool[i] <- sum(p_ord[i, ] < thresh[i, ]) 
  bhglobH0_bool[i] <- any(p_ord[i, ] < thresh[i, ]) 
}                      
mean(bhlokalH0_bool)
mean(bhglobH0_bool)
sum(bhlokalH0_bool)
sum(bhglobH0_bool)
