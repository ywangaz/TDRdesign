##########################################################################
#            Three-regional Dual-criterion One-stage 3 by 2              #
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

getprob3by2 = function(p0,p1,nc,rL,sL,mL,call,LL=1e-7,k=1) {
  # return a probability matrix of shape (nrL,nsL, nmL), where nrL, nsL, and nmL are length of input 
  # parameters vectors under specified response rates under the null and alternative hypothesis
  # rL: vector of integers
  # sL: vector of integers
  # mL: vector of integers
  # nc: size of the control arm
  # k: randomization ratio (E:C), default is 1
  nrL = length(rL)
  nsL = length(sL)
  nmL = length(mL)
  prsm = array(0, c(nrL,nsL, nmL)) # probability matrix of s length rows and m length columns
  ne = k*nc
  
  pec_ = pbinres(p0,p1,nc,k) # given ne and nc, binomial probability of every combination of responses
  
  if ((p1-p0)<0.001) {
    pec = pec_[[2]]
  } else {
    pec = pec_[[1]]
  }
  
  ec1L = which(pec > LL) # reduce computation
  for (ec1 in ec1L) {
    c1 = floor((ec1-0.01)/(ne+1))+1 # number of rows = ne + 1
    e1 = ec1 - (c1-1)*(ne+1)
    z1 = e1 - 1 - k*(c1-1)
    if (call=='power' || call=='alpha') {
      prsm = prsm + outer(outer(999>rL, z1>=sL), (e1-1)>=mL)*pec[e1,c1]
    } else if (call == 'type2' || call=='nogo') {
      prsm = prsm + outer(outer(z1<=rL, 999>sL), 999>mL)*pec[e1,c1] + outer(outer(z1>rL, z1<sL), (e1-1)<mL)*pec[e1,c1]
    } else if (call == 'eta1' || call=='gam1') {
      prsm = prsm + outer(outer(z1>rL, z1<sL), (e1-1)>=mL)*pec[e1,c1]
    } else if (call == 'eta2' || call=='gam2') {
      prsm = prsm + outer(outer(999>rL, z1>=sL), (e1-1)<mL)*pec[e1,c1]
    }
  }
  return(prsm)
}

####---------------------Main Function----------------------####
getTDRfunc1s3by2 = function(x,LL=1e-9,k=1,FILENAME=FILENAME) {
  # x: one entry of experiment data
  # k: randomization ratio (E:C), default is 1
  
  # initiate experiment data
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
  expicmax = x['expicmax'][[1]]
  ci = x['ci'][[1]]

  L=1e-9
  LL=1e-7
  
  # search starting and ending constraints
  # starting sample size use 50% of Hong's sample size
  nc = round(0.5 * ((qnorm(nogomin)*sqrt(2*p0*(1-p0))+qnorm(1-betamax)*sqrt(p0*(1-p0)+p1*(1-p1)))/(p1-p0))^2)
  soln = 0
  ncs = 0
  ncmax = 100
  
  # alpha1 and alpha2 response rates
  p01=p0+.05     
  p11=p1+.05      ;  if(p11>.999) p11=1
  p02=(p0+p1)/2
  p12=p1+(p1-p0)/2  ; if(p12>.999) p12=1
  
  # start searching...
  while (nc <= ncmax) { # nc2 loop step by 4
    if( soln==2 ) { nc = nc+1   }     # first soluton found and no solution found for ne2 
    
    ne = k*nc
    print(ne)
    
    # specify s and m vectors for searching
    se0 = sqrt(p0*(1-p0)*(nc))
    se1 = sqrt(p1*(1-p1)*(ne))
    se10 = sqrt(p1*(1-p1)*(ne)+p0*(1-p0)*(nc))
    
    smin = max(-nc, floor(p0*(ne-k*nc)))
    smax = min(ne, ceiling(p1*ne-p0*k*nc + 4*se10))
    
    rmin = max(-nc, floor(p0*(ne-k*nc)-2*se10))
    rmax = min(ne, ceiling(p1*ne-p0*k*nc + 2*se10))
      
    mmin = floor(p0*(ne))
    mmax = min(ne, ceiling((ne)*p1+4*se1))
    
    rL = rmin:rmax
    sL = smin:smax
    mL = mmin:mmax
    
    p0L = max(0.001, p0-qnorm(0.5+0.5*ci)*sqrt(p0*(1-p0)/((1+k)*(nc))))
    p0U = min(0.999, p0+qnorm(0.5+0.5*ci)*sqrt(p0*(1-p0)/((1+k)*(nc))))
    p0.ci = seq(p0L, p0U, length.out=11)
    alpha.actual = rep(0,11)
    power.actual = rep(0,11)
    
    # control alpha
    alpL = getprob3by2(p0,p0,nc,rL,sL,mL,'alpha',LL,k)
    alpLidx = which(alpL<=alphamax)
    
    # constrain: s < m
    alpLsi = which(alpL <= alphamax, arr.ind = T)
    constrainr = sL[alpLsi[,1]]
    constrains = sL[alpLsi[,2]]
    constrainm = mL[alpLsi[,3]]
    constrain = (constrainr < constrains) & (constrains < constrainm)
    alpLsicontrain = alpLsi[constrain,]
    nr = length(rL)
    ns = length(sL)
    nm = length(mL)
    if (length(alpLsicontrain)==3) {
      constrainx = (alpLsicontrain[3]-1)*nr*ns + (alpLsicontrain[2]-1)*nr+alpLsicontrain[1]
    } else {
    constrainx = (alpLsicontrain[,3]-1)*nr*ns + (alpLsicontrain[,2]-1)*nr+alpLsicontrain[,1]
    }
    xd = intersect(alpLidx, constrainx)
    
    if (length(xd) > 0) { # control beta
      betaL = getprob3by2(p0,p1,nc,rL,sL,mL,'type2',LL,k)
      betaLidx = which(betaL<=betamax)
      xd=intersect(xd, betaLidx)
      
      if (length(xd) > 0) { # control alpha1
        alp1L = getprob3by2(p01,p01,nc,rL,sL,mL,'alpha',LL,k)
        alp1Lidx = which(alp1L <= alpha1max)
        xd = intersect(xd, alp1Lidx)
        
        if (length(xd)>0) { # control alpha2
          alp2L = getprob3by2(p02,p02,nc,rL,sL,mL,'alpha',LL,k)
          alp2Lidx = which(alp2L <= alpha2max)
          xd = intersect(xd, alp2Lidx)
          
          if (length(xd > 0)) { # control alpha max
            for (p0s in p0.ci) {
              alpmaxL = getprob3by2(p0s,p0s,nc,rL,sL,mL,'alpha',LL,k)
              alpmaxLidx = which(alpmaxL<=alphamaxmax)
              xd=intersect(alpmaxLidx,xd)
            }
            if (length(xd) > 0) { # control power
              pwrL = getprob3by2(p0,p1,nc,rL,sL,mL,'power',LL,k)
              pwrLidx = which(pwrL >= pwrmin)
              xd = intersect(pwrLidx, xd)
              
              if (length(xd) > 0) { # control eta
                eta1L = getprob3by2(p0,p0,nc,rL,sL,mL,'eta1',LL,k)
                eta2L = getprob3by2(p0,p0,nc,rL,sL,mL,'eta2',LL,k)
                etaL = eta1L + eta2L
                etaLidx = which(etaL<=etamax)
                xd=intersect(xd, etaLidx)
                
                if (length(xd) > 0) { # control gam
                  gam1L = getprob3by2(p0,p1,nc,rL,sL,mL,'gam1',LL,k)
                  gam2L = getprob3by2(p0,p1,nc,rL,sL,mL,'gam2',LL,k)
                  gamL = gam1L + gam2L
                  gamLidx = which(gamL<=gammax)
                  xd=intersect(xd, gamLidx)
                  
                  if (length(xd) > 0) { # control nogo (accept H0|H0)
                    
                    expicL = (etaL + gamL)/2
                    expicLidx = which(expicL<=expicmax)
                    xd=intersect(expicLidx,xd)
                    
                    nogoL = getprob3by2(p0,p0,nc,rL,sL,mL,'nogo',LL,k)
                    nogoLidx = which(nogoL>=nogomin)
                    xd=intersect(xd, nogoLidx)
                    
                    if (length(xd) > 0) soln = soln + 1
                    } # nogo
                  } # gam
                } # eta
            } # pwr
          } # alphamax
        } # alpha2
      } # alpha1
    } # alpha
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
        midx = floor((iam-0.01)/(nr*ns)) + 1
        sidx = floor((iam - (midx-1)*nr*ns - 0.01)/nr) + 1
        ridx = iam - (midx-1)*nr*ns - (sidx-1)*nr
        r = rL[ridx]
        s = sL[sidx]
        m = mL[midx]
        
        print('solution found')
        alphamaxvalue = -0.01
        for (p0s in p0.ci) { 
          alpmaxL = getprob3by2(p0s,p0s,nc,r,s,m,'alpha',LL,k)
          alphamaxvalue = max(alphamaxvalue, alpmaxL[1])
        }
        
        search_result = c(p0,p1,nc,ne,r,s,m,
                          round(alpL[iam], 5), round(alp1L[iam], 5), round(alp2L[iam], 5),
                          round(alphamaxvalue, 5), round(pwrL[iam], 5),round(betaL[iam], 5),
                          round(etaL[iam], 5),round(gamL[iam], 5), 
                          round(eta1L[iam],5), round(eta2L[iam],5), round(gam1L[iam], 5), round(gam2L[iam], 5), 
                          round(nogoL[iam], 5), round(expicL[iam], 5),
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
tdr_1s_3by2 = function(FILEPREFIX, alphamax, alpha1max, alpha2max, alphamaxmax, betamax, pwrmin, eta, gam, expic, ci=0.9,k=1) {
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
  cat(paste0(c('p0', 'p1', 'nc','ne','r', 's','m',
               'alpha', 'alpha1', 'alpha2', 'alphamax','power','beta',
               'eta', 'gam', 'eta1', 'eta2', 'gam1', 'gam2',
               'nogo','expic', 'N','etamax','gammax','expicmax', 'ci'),collapse = ","))
  closeAllConnections()
  
  # simulations
  TDR_1s_3by2_results = c()
  
  for (i in 1:nrow(dat)) {
    print(i)
    tic()
    s1_3by2 = getTDRfunc1s3by2(dat[i,],k=k, FILENAME=FILENAME)
    print(s1_3by2)
    toc()
    TDR_1s_3by2_results = rbind(TDR_1s_3by2_results, s1_3by2)
  }
  
} 

####-------------------------Parameter Search------------------------####
fileprefix = 'TDR_alpha20_1s_3by2_alphaonly'

for (expic in seq(0.15,0.35,0.05)) {
  for (gam in seq(0.10,0.25,0.01)) {
    tdr_1s_3by2(FILEPREFIX=fileprefix, 
                alphamax=0.20, alpha1max=0.99, alpha2max=0.99,alphamaxmax=0.99, 
                betamax=0.20, pwrmin=0.75,
                eta=0.99, gam=gam, expic=expic, ci=0.90,
                k=1)
  }
}

fileprefix = 'TDR_alpha10_1s_3by2_alphaonly'

for (expic in seq(0.15,0.35,0.05)) {
  for (gam in seq(0.10,0.25,0.01)) {
    tdr_1s_3by2(FILEPREFIX=fileprefix, 
                alphamax=0.10, alpha1max=0.99, alpha2max=0.99,alphamaxmax=0.99, 
                betamax=0.10, pwrmin=0.85,
                eta=0.99, gam=gam, expic=expic, ci=0.90,
                k=1)
  }
}

fileprefix = 'TDR_alpha20_1s_3by2_alphamax30'

for (expic in seq(0.15,0.35,0.05)) {
  for (gam in seq(0.10,0.25,0.01)) {
    tdr_1s_3by2(FILEPREFIX=fileprefix, 
                alphamax=0.20, alpha1max=0.99, alpha2max=0.99,alphamaxmax=0.20, 
                betamax=0.20, pwrmin=0.75,
                eta=0.99, gam=gam, expic=expic, ci=0.30,
                k=1)
  }
}

fileprefix = 'TDR_alpha10_1s_3by2_alphamax30'

for (expic in seq(0.15,0.35,0.05)) {
  for (gam in seq(0.10,0.25,0.01)) {
    tdr_1s_3by2(FILEPREFIX=fileprefix, 
                alphamax=0.10, alpha1max=0.99, alpha2max=0.99,alphamaxmax=0.10, 
                betamax=0.10, pwrmin=0.85,
                eta=0.99, gam=gam, expic=expic, ci=0.30,
                k=1)
  }
}

##########################################################################
#                                  End                                   #
##########################################################################

