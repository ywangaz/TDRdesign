##########################################################################
#            Three-regional Dual-criterion One-stage 2 by 2              #
##########################################################################
# -----------------------------------------------------------------------#
# code reference: Litwin et al. 2007
# -----------------------------------------------------------------------#

library(dplyr)
library(tictoc)

pbinres <- function(p0, p1, nc, k=1) {
  # return probability matrix of possible combinations of responses in the experiemental arm and
  # the control arm under the null and alternative hypothesis
  # nc: size of the control arm
  # k: randomization ratio (E:C), default is 1
  ne=nc*k
  ema=matrix( dbinom(0:ne,ne,p1), ne+1, 1 ) # +1 because of #res starts from 0 >> n+1 possibility
  emn=matrix( dbinom(0:ne,ne,p0), ne+1, 1 )
  cm=matrix( dbinom(0:nc,nc,p0), 1, nc+1 )
  ecma=ema%*%cm # matrix of all possible combination of responses under alternative hypothesis
  ecmn=emn%*%cm # matrix of all possible combination of responses under null hypothesis
  return( list(ecma,ecmn ) )      #  note if x=pbinres(p0,p1,ne) then ecma=x[[1]]  ecmn=x[[2]]
}

getprob2by2 <- function(p0, p1, nc, sL, mL, call, L=1e-7, k=1) {
  # return a probability matrix of shape (nsL, nmL), where nsL and nmL are length of input 
  # parameters vectors under specified response rates under the null and alternative hypothesis
  # sL: vector of integers
  # mL: vector of integers
  # nc: size of the control arm
  # k: randomization ratio (E:C), default is 1
  nsL = length(sL)
  nmL = length(mL)
  psm = array(0, c(nsL, nmL)) # probability matrix of s length rows and m length columns
  ne = k*nc
  
  pec_ = pbinres(p0,p1,nc,k) # given ne and nc, binomial probability of every combination of responses
  
  if ((p1-p0)<0.001) {
    pec = pec_[[2]]
  } else {
    pec = pec_[[1]]
  }
  
  ec1L = which(pec > L) # reduce computation
  for (ec1 in ec1L) {
    c1 = floor((ec1-0.01)/(ne+1))+1 # number of rows = ne + 1
    e1 = ec1 - (c1-1)*(ne+1)
    z1 = e1 - 1 - k*(c1-1)
    if (call=='power' || call=='alpha') {
      psm = psm + outer(z1>=sL, (e1-1)>=mL)*pec[e1,c1]
    } else if (call == 'type2' || call=='nogo') {
      psm = psm + outer(z1<sL, 999>mL)*pec[e1,c1]
    } else if (call == 'eta' || call=='gam') {
      psm = psm + outer(z1>=sL, (e1-1)<mL)*pec[e1,c1]
    }
  }
  return(psm)
}

####---------------------Main Function----------------------####
getTDRfunc1s2by2 <- function(x,LL=1e-9,k=1,FILENAME=FILENAME) {
  # x: one entry of experiment data
  # k: randomization ratio (E:C), default is 1
  
  # initiate experiment data
  p0=x['p0'][[1]]
  p1=x['p1'][[1]]
  alphamax = x['alphamax'][[1]]
  alpha1max = x['alpha1max'][[1]]
  alpha2max = x['alpha2max'][[1]]
  alphamaxmax = x['alphamaxmax'][[1]]
  ci = x['ci'][[1]]
  betamax = x['betamax'][[1]]
  pwrmin = x['pwrmin'][[1]]
  nogomin = x['nogomin'][[1]]
  etamax = x['etamax'][[1]]
  gammax = x['gammax'][[1]]
  expicmax = x['expicmax'][[1]]
  
  # search starting and ending constraints
  nc = 10
  ncs = 0
  ncmax = 150
  soln = 0
  
  # alpha1 and alpha2 response rates
  p01=p0+.05     
  p11=p1+.05      ;  if(p11>.999) p11=1
  p02=(p0+p1)/2
  p12=p1+(p1-p0)/2  ; if(p12>.999) p12=1
  
  # start searching...
  while (nc<=ncmax) {
    print(soln)
    if( soln==2 ) { nc = nc+1   } 
    ne = k*nc
    print(ne)
    
    # specify s and m vectors for searching
    se1 = sqrt(p1*(1-p1)*(ne))
    se10 = sqrt(p1*(1-p1)*(ne)+p0*(1-p0)*(nc))
    
    smin = max(-nc, floor(p0*(ne-k*nc)))
    smax = min(ne, ceiling(p1*ne-p0*k*nc + 4*se10))
    
    mmin = floor(p0*(ne))
    mmax = min(ne, ceiling((ne)*p1+4*se1))
    
    sL = smin:smax
    mL = mmin:mmax
    
    p0L = max(0.001, p0-qnorm(0.5+0.5*ci)*sqrt(p0*(1-p0)/((1+k)*(nc))))
    p0U = min(0.999, p0+qnorm(0.5+0.5*ci)*sqrt(p0*(1-p0)/((1+k)*(nc))))
    p0.ci = seq(p0L, p0U, length.out=11)
    
    # control alpha
    alpL = getprob2by2(p0,p0,nc,sL,mL,'alpha')
    alpLidx = which(alpL<=alphamax)
    
    # constrain: s < m
    alpLsi = which(alpL <= alphamax, arr.ind = T)
    constrains = sL[alpLsi[,1]]
    constrainm = mL[alpLsi[,2]]
    constrain = constrains < constrainm
    alpLsicontrain = alpLsi[constrain,]
    ns = length(sL)
    nm = length(mL)
    if (length(alpLsicontrain)==2) {
      constrainx = (alpLsicontrain[2]-1)*ns + alpLsicontrain[1]
    } else {
      constrainx = (alpLsicontrain[,2]-1)*ns + alpLsicontrain[,1]
    }
    xd = intersect(alpLidx, constrainx)
    
    if (length(xd) > 0) {     # control beta
      betaL = getprob2by2(p0,p1,nc,sL,mL,'type2')
      betaLidx = which(betaL<=betamax)
      xd=intersect(xd, betaLidx)
      
      if (length(xd) > 0) { # control alpha1
        alp1L = getprob2by2(p01,p01,nc,sL,mL,'alpha')
        alp1Lidx = which(alp1L <= alpha1max)
        xd = intersect(xd, alp1Lidx)
        
        if (length(xd)>0) { # control alpha2
          alp2L = getprob2by2(p02,p02,nc,sL,mL,'alpha')
          alp2Lidx = which(alp2L <= alpha2max)
          xd = intersect(xd, alp2Lidx)
          
          if (length(xd > 0)) { # control alpha max
            for (p0s in p0.ci) {
              alpmaxL = getprob2by2(p0s,p0s,nc,sL,mL,'alpha')
              alpmaxLidx = which(alpmaxL<=alphamaxmax)
              xd=intersect(alpmaxLidx,xd)
            }
            if (length(xd) > 0) { # control power
              pwrL = getprob2by2(p0,p1,nc,sL,mL,'power')
              pwrLidx = which(pwrL >= pwrmin)
              xd=intersect(xd, pwrLidx)
              
              if (length(xd) > 0) { # control eta
                etaL = getprob2by2(p0,p0,nc,sL,mL,'eta')
                etaLidx = which(etaL<=etamax)
                xd=intersect(xd, etaLidx)
                
                if (length(xd) > 0) { # control gam
                  gamL = getprob2by2(p0,p1,nc,sL,mL,'gam')
                  gamLidx = which(gamL<=gammax)
                  xd=intersect(xd, gamLidx)

                  expicL = (etaL + gamL)/2
                  expicLidx = which(expicL<=expicmax)
                  xd=intersect(expicLidx,xd)
                  
                  nogoL = getprob2by2(p0,p0,nc,sL,mL,'nogo')
                  nogoLidx = which(nogoL>=nogomin)
                  xd=intersect(xd, nogoLidx)
                  
                  if (length(xd) > 0) soln = soln + 1
                } # gam, expic, nogo
              } # eta
            } # power
          } # control alphamax
        } # control alpha2
      } # control alpha1
    } # control beta
    
    if( soln==0 ) { nc=nc+4 }  
    
    if (soln == 1) {
      ncs = nc
      nc = max(1,nc-5)
      soln = 2
    }
    
    if (soln==3 || nc==ncs) {
      
      jz = which(alpL == min(alpL[xd]))
      jz = intersect(xd, jz)
      iam = jz[1]
      if (length(iam) < 0) {
        iam = xd[1]
      } else {
          midx = floor((iam-0.01)/ns) + 1
          sidx = iam - (midx-1)*ns
          s = sL[sidx]
          m = mL[midx]
          
          print('solution found')
          
          alphamaxvalue = -0.01
          for (p0s in p0.ci) { 
            alpmaxL = getprob2by2(p0s,p0s,nc,s,m,'alpha')
            alphamaxvalue = max(alphamaxvalue, alpmaxL[1])
          }
          
          search_result = c(p0,p1,nc,ne,s,m,
                            round(alpL[iam], 5), round(alp1L[iam], 5), round(alp2L[iam], 5),round(alphamaxvalue, 5), 
                            round(pwrL[iam], 5),round(betaL[iam], 5),
                            round(etaL[iam], 5),round(gamL[iam], 5),round(nogoL[iam], 5), round(expicL[iam], 5),
                            nc+ne,etamax, gammax, expicmax, ci)
          
          sink(paste('results/',FILENAME, sep=''), append=T)
          cat('\n')
          cat(paste0(search_result, collapse = ","))
          closeAllConnections()
          
          return(search_result)
          break
        }
    }
  } # while loop
}

####---------------------Experiment---------------------####
tdr_1s_2by2 <- function(FILEPREFIX, alphamax, alpha1max, alpha2max, alphamaxmax, betamax, pwrmin, eta, gam, expic, ci=0.9,k=1) {
  FILENAME = paste(FILEPREFIX,'_eta', round(eta*100),'_gam', round(gam*100), '_ic',round(expic*100), '_confidence',round(ci*100), '.csv', 
                   sep='')
  print(FILENAME)
  etamax = eta
  gammax = gam
  expicmax = expic
  ci = ci
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
                   etamax=etamax,
                   gammax=gammax,
                   expicmax=expicmax,
                   ci=ci)
  
  # write data
  sink(paste('results/',FILENAME, sep=''), append=F)
  cat('\n')
  cat(FILENAME)
  cat('\n')
  cat(paste0(c('p0', 'p1', 'nc','ne','s','m',
               'alpha', 'alpha1', 'alpha2', 'alphamax','power','beta','eta', 'gam', 'nogo','expic', 'N',
               'etamax','gammax','expicmax', 'ci'),collapse = ","))
  closeAllConnections()
  
  # run simulations
  TDR_1s_2by2_results = c()
  for (i in 1:nrow(dat)) {
    print(i)
    tic()
    s1 = getTDRfunc1s2by2(dat[i,],k=k, FILENAME=FILENAME)
    print(s1)
    toc()
    TDR_1s_2by2_results = rbind(TDR_1s_2by2_results, s1)
  }
  
}

####-------------------------Parameter Search------------------------####
  
fileprefix = 'TDR_alpha20_1s_2by2_alphaonly'

for (expic in seq(0.05,0.30,0.05)) {
  for (gam in seq(0.01,0.2,0.01)) {
    tdr_1s_2by2(FILEPREFIX=fileprefix, 
                alphamax=0.20, alpha1max=0.99, alpha2max=0.99,alphamaxmax=0.99, 
                betamax=0.20, pwrmin=0.75,
                eta=0.99, gam=gam, expic=expic, ci=0.90,
                k=1)
  }
}

fileprefix = 'TDR_alpha10_1s_2by2_alphaonly'

for (expic in seq(0.05,0.30,0.05)) {
  for (gam in seq(0.01,0.2,0.01)) {
    tdr_1s_2by2(FILEPREFIX=fileprefix, 
                alphamax=0.10, alpha1max=0.99, alpha2max=0.99,alphamaxmax=0.99, 
                betamax=0.10, pwrmin=0.85,
                eta=0.99, gam=gam, expic=expic, ci=0.90,
                k=1)
  }
}

fileprefix = 'TDR_alpha20_1s_2by2_alphamax30'

for (expic in seq(0.05,0.30,0.05)) {
  for (gam in seq(0.01,0.2,0.01)) {
    tdr_1s_2by2(FILEPREFIX=fileprefix, 
                alphamax=0.20, alpha1max=0.99, alpha2max=0.99,alphamaxmax=0.20, 
                betamax=0.20, pwrmin=0.75,
                eta=0.99, gam=gam, expic=expic, ci=0.30,
                k=1)
  }
}

fileprefix = 'TDR_alpha10_1s_2by2_alphamax30'

for (expic in seq(0.05,0.30,0.05)) {
  for (gam in seq(0.01,0.2,0.01)) {
    tdr_1s_2by2(FILEPREFIX=fileprefix, 
                alphamax=0.10, alpha1max=0.99, alpha2max=0.99,alphamaxmax=0.10, 
                betamax=0.10, pwrmin=0.85,
                eta=0.99, gam=gam, expic=expic, ci=0.30,
                k=1)
  }
}


##########################################################################
#                                  End                                   #
##########################################################################