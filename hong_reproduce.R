#### Hong & Wang Design ####
# code reference: Litwin et. al. 2007

library(dplyr)
library(tictoc)
library(ph2bayes)

pbinreshong=function(pc,pe,nc,k=1) {
  ne=nc*k # 1:1 randomization
  p = pc #(pc+pe)/2
  
  ema=matrix(dbinom(0:ne,ne,pe), ne+1, 1 ) # +1 because of #res starts from 0 >> n+1 possibility
  cm=matrix(dbinom(0:nc,nc,pc), 1, nc+1 )
  
  emn2 = matrix(dbinom(0:ne,ne,p), ne+1,1)
  cm2 = matrix(dbinom(0:nc,nc,p), 1, nc+1 )
  ecma=ema%*%cm # matrix of all possible combination of responses under alternative hypothesis
  ecmn=emn2%*%cm2 # matrix of all possible combination of responses under null hypothesis
  return( list(ecma,ecmn ) )      #  note if x=pbinreshong(pc,pe,ne) then ecma=x[[1]]  ecmn=x[[2]]
  
}

getprobHong = function(p0,p1,nc,d1L,d2L,call,LL=1e-7,k=1) {
  
  nd1L = length(d1L)
  nd2L = length(d2L)
  pd1d2 = array(0, c(nd1L, nd2L)) # probability matrix of s length rows and m length columns
  ne = k*nc
  
  pec_ = pbinreshong(p0,p1,nc,k) # given ne and nc, binomial probability of every combination of responses
  
  if ((p1-p0)<0.001) {
    pec = pec_[[2]]
  } else {
    pec = pec_[[1]]
  }
  
  # ec1L = which(pec > 0) # reduce computation
  for (e in 0:ne) { 
    for(c in 0:nc) {
      z1 = e/ne - c/nc # observed response rate difference
      if (call=='power' || call=='alpha') {
        pd1d2 = pd1d2 + outer(999>=d1L, z1>=d2L)*pec[e+1,c+1]
      } else if (call == 'type2' || call=='nogo') {
        pd1d2 = pd1d2 + outer(z1<=d1L, 999>d2L)*pec[e+1,c+1]
      } else if (call == 'eta' || call=='gam') {
        pd1d2 = pd1d2 + outer(z1>d1L, z1<d2L)*pec[e+1,c+1]
      }
    }
  }
  return(pd1d2)
}

# getprobHong(0.10,0.30,51,0.0778,0.137,'type2')
# getprobHong(0.10,0.30,51,0.095,0.1365,'gam')
# getprobHong(0.20,0.20,51,0.078,0.137,'alpha')
# getprobHong(0.10,0.10,51,0.078,0.137,'eta')

getTDRfunc1sHong = function(x,LL=1e-9,k=1,FILENAME=FILENAME) {
  
  p0=x['p0'][[1]]
  p1=x['p1'][[1]]
  alphamax = x['alphamax'][[1]]
  alpha1max = x['alpha1max'][[1]]
  alpha2max = x['alpha2max'][[1]]
  pwrmin = x['pwrmin'][[1]]
  alphamaxmax = x['alphamaxmax'][[1]]
  betamax = x['betamax'][[1]]
  nogomin = x['nogomin'][[1]]
  etamax = x['etamax'][[1]]
  gammax = x['gammax'][[1]]
  ci = x['ci'][[1]]
  
  soln = 0
  
  nc = 10
  ncs = 0
  ncmax = 200
  soln = 0
  
  p01=p0+.05     
  p11=p1+.05      ;  if(p11>.999) p11=1
  p02=(p0+p1)/2
  p12=p1+(p1-p0)/2  ; if(p12>.999) p12=1

  while (nc<=ncmax) {
    print(soln)
    if( soln==2 ) { nc = nc+1   } 
    ne = k*nc
    print(ne)
  
    d1L = seq(0.05,0.10,0.001)
    d2L = seq(0.09,0.20,0.001)
    
    p0L = max(0.001, p0-qnorm(0.5+0.5*ci)*sqrt(p0*(1-p0)/((1+k)*(nc))))
    p0U = min(0.999, p0+qnorm(0.5+0.5*ci)*sqrt(p0*(1-p0)/((1+k)*(nc))))
    p0.ci = seq(p0L, p0U, length.out=11)
    alpha.actual = rep(0,11)
    
    # control alpha
    alpL = getprobHong(p0,p0,nc,d1L,d2L,'alpha')
    alpLidx = which(alpL<=alphamax)
    
    # constrain: s < m
    alpLsi = which(alpL <= alphamax, arr.ind = T)
    constraind1 = d1L[alpLsi[,1]]
    constraind2 = d2L[alpLsi[,2]]
    constrain = constraind1 < constraind2
    alpLsicontrain = alpLsi[constrain,]
    nd1 = length(d1L)
    nd2 = length(d2L)
    if (length(alpLsicontrain)==2) {
      constrainx = (alpLsicontrain[2]-1)*nd1 + alpLsicontrain[1]
    } else {
      constrainx = (alpLsicontrain[,2]-1)*nd1 + alpLsicontrain[,1]
    }
    xd = intersect(alpLidx, constrainx)
    
    # control beta
    if (length(xd) > 0) {
      betaL = getprobHong(p0,p1,nc,d1L,d2L,'type2')
      betaLidx = which(betaL<=betamax)
      xd = intersect(xd, betaLidx)
      
      if (length(xd) > 0) { # control alpha1
        alp1L = getprobHong(p01,p01,nc,d1L,d2L,'alpha')
        alp1Lidx = which(alp1L <= alpha1max)
        xd = intersect(xd, alp1Lidx)
        
        if (length(xd)>0) { # control alpha2
          alp2L = getprobHong(p02,p02,nc,d1L,d2L,'alpha')
          alp2Lidx = which(alp2L <= alpha2max)
          xd = intersect(xd, alp2Lidx)
          
          if (length(xd > 0)) {
            # print(c(max(pwrL[xd]), min(alpL[xd])))
            # control alpha max
            for (p0s in p0.ci) { # control alpha max
              alpmaxL = getprobHong(p0s,p0s,nc,d1L,d2L,'alpha')
              alpmaxLidx = which(alpmaxL<=alphamaxmax)
              xd=intersect(alpmaxLidx,xd)
            }
            
            if (length(xd > 0)) { # control power
              
              # control power
              pwrL = getprobHong(p0,p1,nc,d1L,d2L,'power')
              pwrLidx = which(pwrL >= pwrmin)
              xd=intersect(xd, betaLidx)
              
              if (length(xd) > 0) { # control eta
                etaL = getprobHong(p0,p0,nc,d1L,d2L,'eta')
                etaLidx = which(etaL<=etamax)
                xd=intersect(xd, etaLidx)
                # print(c(max(pwrL[xd]), min(alpL[xd]), min(betaL[xd])))
                
                if (length(xd) > 0) { # control nogo (accept H0|H0)
                  gamL = getprobHong(p0,p1,nc,d1L,d2L,'gam')
                  gamLidx = which(gamL<=gammax)
                  xd=intersect(xd, gamLidx)
                  
                  nogoL = getprobHong(p0,p0,nc,d1L,d2L,'nogo')
                  # print(c(max(pwrL[xd]), min(alpL[xd]), min(betaL[xd]),max(nogoL[xd])))
                  
                  if (length(xd) > 0) soln = soln + 1
                } # gamma
              } # eta
            } # beta
          } # control alphamax
        } # control alpha2
      } # control alpha1
    } # alpha
    
    if( soln==0 ) { nc=nc+4 }  
    
    if (soln == 1) {
      ncs = nc
      nc = max(1,nc-5)
      soln = 2
    }
    
    if (soln==3 || nc==ncs) {
      jz = which(alpL==min(alpL[xd])) 
      jz = intersect(xd, jz)
      iam = jz[1]
      if (length(iam) <= 0) {
        iam = xd[1]
      }
      if (length(iam) > 0) {
        d1idx = floor((iam-0.01)/nd2) + 1
        d2idx = iam - (d1idx-1)*nd2
        d1 = d1L[d1idx]
        d2 = d2L[d2idx]
        print('solution found')
        
        alphamaxvalue = -0.01
        for (p0s in p0.ci) { 
          alpmaxL = alpmaxL = getprobHong(p0s,p0s,nc,d1,d2,'alpha')
          alphamaxvalue = max(alphamaxvalue, alpmaxL[1])
        }

        search_result = c(p0,p1,nc,ne,d1,d2,
                          round(alpL[iam], 5), round(alp1L[iam], 5), round(alp2L[iam], 5),
                          round(alphamaxvalue, 5), round(pwrL[iam], 5),round(betaL[iam], 5),
                          round(etaL[iam], 5),round(gamL[iam], 5),round(nogoL[iam], 5), nc+ne)
        print(search_result)
        
        sink(paste('results/',FILENAME, sep=''), append=T)
        cat('\n')
        cat(paste0(search_result, collapse = ","))
        closeAllConnections()
        
        return(search_result)
        break
      }
    } # hong soln found
    
  } # while loop
}

#### simulations ####

hong_1s = function(FILEPREFIX, alphamax, alpha1max, alpha2max, alphamaxmax, betamax, pwrmin, eta, gam, ci=0.9,k=1) {
  FILENAME = paste(FILEPREFIX,'_eta', round(eta*100),'_gam', round(gam*100), '_confidence',round(ci*100), '.csv', 
                   sep='')
  print(FILENAME)
  nogomin = 0.01
  
  dat = data.frame(p0=rep(seq(0.05,0.70,0.05),4),
                   p1=c(seq(0.15,0.80,0.05), seq(0.20,0.85,0.05),seq(0.25,0.90,0.05),seq(0.30,0.95,0.05)),
                   es0min=0.50,
                   es1max=0.05,
                   alphamax=alphamax,
                   betamax=betamax,
                   alpha1max=alpha1max,
                   alpha2max=alpha2max,
                   alphamaxmax=alphamaxmax,
                   pwrmin=pwrmin,
                   nogomin=nogomin,
                   etamax=eta,
                   gammax=gam,
                   # expicmax=expicmax,
                   ci=ci)
  
  # write data
  sink(paste('results/',FILENAME, sep=''), append=F)
  cat('\n')
  cat(FILENAME)
  cat('\n')
  cat(paste0(c('p0', 'p1', 'nc','ne','d1','d2',
               'alpha', 'alpha1', 'alpha2', 'alphamax','power','beta','eta', 'gam', 'nogo','N'),collapse = ","))
  closeAllConnections()
  
  TDR_1s_hong_results = c()
  for (i in 1:nrow(dat)) {
    print(i)
    tic()
    hong1 = getTDRfunc1sHong(dat[i,],k=k, FILENAME=FILENAME)
    print(hong1)
    toc()
    TDR_1s_hong_results = rbind(TDR_1s_hong_results, hong1)
  }
  
  # colnames(TDR_1s_hong_results) = c('p0', 'p1', 'nc','ne','d1','d2',
  #                                   'alpha', 'alpha1', 'alpha2', 'alphamax','power','beta','eta', 'gam', 'nogo','N')
  # rownames(TDR_1s_hong_results) = 1:nrow(TDR_1s_hong_results)
  
}

#### experiment ####
FILENAME = ''

# fileprefix = 'hong_alpha20_alphaonly_pwr75'
# hong_1s(FILEPREFIX=fileprefix, 
#         alphamax=0.20, alpha1max=0.99, alpha2max=0.99,alphamaxmax=0.99, 
#         betamax=0.20, pwrmin=0.75,
#         eta=0.10, gam=0.10, ci=0.90,
#         k=1)
# 
# fileprefix = 'hong_alpha20_alphaonly_pwr75'
# hong_1s(FILEPREFIX=fileprefix, 
#         alphamax=0.20, alpha1max=0.99, alpha2max=0.99,alphamaxmax=0.99, 
#         betamax=0.20, pwrmin=0.75,
#         eta=0.20, gam=0.20, ci=0.90,
#         k=1)
# 
# fileprefix = 'hong_alpha10_alphaonly_pwr85'
# hong_1s(FILEPREFIX=fileprefix, 
#         alphamax=0.10, alpha1max=0.99, alpha2max=0.99,alphamaxmax=0.99, 
#         betamax=0.10, pwrmin=0.85,
#         eta=0.10, gam=0.10, ci=0.90,
#         k=1)

fileprefix = 'hong_alpha20_alphamax30'
hong_1s(FILEPREFIX=fileprefix, 
        alphamax=0.20, alpha1max=0.99, alpha2max=0.99,alphamaxmax=0.20, 
        betamax=0.20, pwrmin=0.00,
        eta=0.10, gam=0.10, ci=0.30,
        k=1)

# fileprefix = 'hong_alpha20_alphamax30_pwr75'
# hong_1s(FILEPREFIX=fileprefix, 
#         alphamax=0.20, alpha1max=0.99, alpha2max=0.99,alphamaxmax=0.20, 
#         betamax=0.20, pwrmin=0.75,
#         eta=0.10, gam=0.10, ci=0.30,
#         k=1)

# fileprefix = 'hong_alpha20_alphamax30_pwr75'
# hong_1s(FILEPREFIX=fileprefix, 
#         alphamax=0.20, alpha1max=0.99, alpha2max=0.99,alphamaxmax=0.20, 
#         betamax=0.20, pwrmin=0.75,
#         eta=0.20, gam=0.20, ci=0.30,
#         k=1)

# fileprefix = 'hong_alpha10_alphamax30_pwr85'
# hong_1s(FILEPREFIX=fileprefix, 
#         alphamax=0.10, alpha1max=0.99, alpha2max=0.99,alphamaxmax=0.10, 
#         betamax=0.10, pwrmin=0.85,
#         eta=0.10, gam=0.10, ci=0.30,
#         k=1)





















