

dat = data.frame(p0=rep(seq(0.05,0.70,0.05),4),
                 p1=c(seq(0.15,0.80,0.05), seq(0.20,0.85,0.05),seq(0.25,0.90,0.05),seq(0.30,0.95,0.05)))
dat$p01 = dat$p0 + 0.05
dat$p02 = (dat$p0 + dat$p1)/2

# 1:1 randomization
alpha = 0.20
beta = 0.20
n_alphaonly = 2*ceiling(((qnorm(1-alpha)+qnorm(1-beta))/(dat$p1 - dat$p0))^2*(dat$p0*(1-dat$p0)+dat$p1*(1-dat$p1)))
n_alpha1 = 2*ceiling(((qnorm(1-alpha)+qnorm(1-beta))/(dat$p1 - dat$p01))^2*(dat$p01*(1-dat$p01)+dat$p1*(1-dat$p1)))
n_alpha2 = 2*ceiling(((qnorm(1-alpha)+qnorm(1-beta))/(dat$p1 - dat$p02))^2*(dat$p02*(1-dat$p02)+dat$p1*(1-dat$p1)))
conventional = cbind(dat,n_alphaonly,n_alpha1,n_alpha2)
write.csv(conventional, 'results/conventional_alpha20_1to1.csv')

alpha = 0.10
beta = 0.10
n_alphaonly = 2*ceiling(((qnorm(1-alpha)+qnorm(1-beta))/(dat$p1 - dat$p0))^2*(dat$p0*(1-dat$p0)+dat$p1*(1-dat$p1)))
n_alpha1 = 2*ceiling(((qnorm(1-alpha)+qnorm(1-beta))/(dat$p1 - dat$p01))^2*(dat$p01*(1-dat$p01)+dat$p1*(1-dat$p1)))
n_alpha2 = 2*ceiling(((qnorm(1-alpha)+qnorm(1-beta))/(dat$p1 - dat$p02))^2*(dat$p02*(1-dat$p02)+dat$p1*(1-dat$p1)))
conventional = cbind(dat,n_alphaonly,n_alpha1,n_alpha2)
write.csv(conventional, 'results/conventional_alpha10_1to1.csv')

# n_alphaonly = ceiling((qnorm(1-alpha)*sqrt(2*dat$p0*(1-dat$p0))+qnorm(1-beta)*sqrt(dat$p1*(1-dat$p1)+dat$p0*(1-dat$p0)))^2/((dat$p1-dat$p0)^2))*2
# dat$p01 = dat$p0 + 0.05
# dat$p02 = (dat$p0 + dat$p1)/2
# n_alpha1 = ceiling((qnorm(1-alpha)*sqrt(2*dat$p01*(1-dat$p01))+qnorm(1-beta)*sqrt(dat$p1*(1-dat$p1)+dat$p01*(1-dat$p01)))^2/((dat$p1-dat$p01)^2))*2
# n_alpha2 = ceiling((qnorm(1-alpha)*sqrt(2*dat$p02*(1-dat$p02))+qnorm(1-beta)*sqrt(dat$p1*(1-dat$p1)+dat$p02*(1-dat$p02)))^2/((dat$p1-dat$p02)^2))*2
# conventional = cbind(dat,n_alphaonly,n_alpha1,n_alpha2)
