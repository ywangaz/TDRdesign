library(readr)
hdat <- read_csv("results/hong_control_alphaonly.csv", 
                                   skip = 2)

n_alphaonly=c()
for (i in 1:nrow(hdat)) {
  a = hdat$alpha[i]
  
  pwr = hdat$power[i]
  n_alphaonly[i] = ceiling((qnorm(1-a)*sqrt(2*hdat$p0[i]*(1-hdat$p0[i]))+qnorm(pwr)*sqrt(hdat$p1[i]*(1-hdat$p1[i])+hdat$p0[i]*(1-hdat$p0[i])))^2/((hdat$p1[i]-hdat$p0[i])^2))
}
hdat = cbind(hdat, n_alphaonly)
hdat$n_alphaonly_diff = hdat$n_alphaonly - hdat$ne
write.csv(hdat, 'results/hongvsconvention_alphaonly.csv')
