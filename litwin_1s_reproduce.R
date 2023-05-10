#### Litwin 1s 2by2 ####
# code reference: Litwin et. al. 2007

library(dplyr)
library(tictoc)
library(ph2bayes)

ecmf=function(pc,pe,nc,k=1) {
  ne=nc*k # 1:1 randomization
  ema=matrix( dbinom(0:ne,ne,pe), ne+1, 1 ) # +1 because of #res starts from 0 >> n+1 possibility
  emn=matrix( dbinom(0:ne,ne,pc), ne+1, 1 ) 
  cm=matrix( dbinom(0:nc,nc,pc), 1, nc+1 )
  ecma=ema%*%cm # matrix of all possible combination of responses under alternative hypothesis
  ecmn=emn%*%cm # matrix of all possible combination of responses under null hypothesis
  return( list(ecma,ecmn ) )      #  note if x=ecmf(pc,pe,ne) then ecma=x[[1]]  ecmn=x[[2]]
  
}

getlitwinprob2by2 = function(p0,p1,nc,sL,mL,call,LL=1e-7,k=1) {
  nsL = length(sL)
  nmL = length(mL)
  psm = array(0, c(nsL, nmL)) # probability matrix of s length rows and m length columns
  ne = k*nc
  
  pec_ = ecmf(p0,p1,nc,k) # given ne and nc, binomial probability of every combination of responses
  
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
      psm = psm + outer(z1>=sL, (e1-1)>=mL)*pec[e1,c1]
    } else if (call == 'type2' || call=='nogo') {
      psm = psm + outer(z1<sL, 999>mL)*pec[e1,c1] + outer(z1>=sL, (e1-1)<mL)*pec[e1,c1] # litwin type 2 error
    } else if (call == 'eta' || call=='gam') {
      psm = psm + outer(z1>=sL, (e1-1)<mL)*pec[e1,c1]
    }
  }
  return(psm)
}

getlitwinfunc1s = function(x,LL=1e-9,k=1, FILENAME=FILENAME) {
  
  p0=x['p0'][[1]]
  p1=x['p1'][[1]]
  alphamax = x['alphamax'][[1]]
  alpha1max = x['alpha1max'][[1]]
  alpha2max = x['alpha2max'][[1]]
  pwrmin = x['pwrmin'][[1]]
  alphamaxmax = x['alphamaxmax'][[1]]
  betamax = x['betamax'][[1]]
  nogomin = x['nogomin'][[1]]
  # expicmax = x['expicmax'][[1]]
  ci = x['ci'][[1]]
  nc = 10
  ncs = 0
  ncmax = 150
  soln = 0
  
  p01=p0+.05     
  p11=p1+.05      ;  if(p11>.999) p11=1
  p02=(p0+p1)/2
  p12=p1+(p1-p0)/2  ; if(p12>.999) p12=1
  
  while (nc<=ncmax) {
    # print(soln)
    if( soln==2 ) { nc = nc+1   } 
    ne = k*nc
    print(ne)
    
    se0 = sqrt(p0*(1-p0)*(nc))
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
    alpha.actual = rep(0,11)
    power.actual = rep(0,11)
    
    # control alpha
    alpL = getlitwinprob2by2(p0,p0,nc,sL,mL,'alpha')
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
    
    # control beta
    if (length(xd) > 0) {
      betaL = getlitwinprob2by2(p0,p1,nc,sL,mL,'type2')
      betaLidx = which(betaL<=betamax)
      xd = intersect(xd, betaLidx)
      
      if (length(xd) > 0) { # control alpha1
        alp1L = getlitwinprob2by2(p01,p01,nc,sL,mL,'alpha')
        alp1Lidx = which(alp1L <= alpha1max)
        xd = intersect(xd, alp1Lidx)
        
        if (length(xd)>0) {
          alp2L = getlitwinprob2by2(p02,p02,nc,sL,mL,'alpha')
          alp2Lidx = which(alp2L <= alpha2max)
          xd = intersect(xd, alp2Lidx)
          
          if (length(xd > 0)) {
            for (p0s in p0.ci) { # control alpha max
              alpmaxL = getlitwinprob2by2(p0s,p0s,nc,sL,mL,'alpha')
              alpmaxLidx = which(alpmaxL<=alphamaxmax)
              xd=intersect(alpmaxLidx,xd)
            }
            if (length(xd) > 0) { # control power
              pwrL = getlitwinprob2by2(p0,p1,nc,sL,mL,'power')
              pwrLidx = which(pwrL >= pwrmin)
              xd=intersect(xd, pwrLidx)
              # print(c(max(pwrL[xd]), min(alpL[xd]), min(betaL[xd])))
              
              if (length(xd) > 0) { # control nogo (accept H0|H0)
                nogoL = getlitwinprob2by2(p0,p0,nc,sL,mL,'nogo')
                nogoLidx = which(nogoL>=nogomin)
                xd=intersect(xd, nogoLidx)
                # print(c(max(pwrL[xd]), min(alpL[xd]), min(betaL[xd]),max(nogoL[xd])))
                
                if (length(xd) > 0) soln = soln + 1
              } # nog0
            } # power
          } # alpha max
        } # alpha2
      } # alpha1
    } # beta
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
            alpmaxL = getlitwinprob2by2(p0s,p0s,nc,s,m,'alpha')
            alphamaxvalue = max(alphamaxvalue, alpmaxL[1])
          }
          
          eta = getlitwinprob2by2(p0,p0,nc,s,m,'eta')
          gam = getlitwinprob2by2(p0,p1,nc,s,m,'gam')
          # calculated expected probability of inconclusive
          expic = 0
          p0e_ = seq(p0,p1,0.01) # assume uniform distribution between p0 and p1
          np0e = length(p0e_)
          for (p0e in p0e_) {
            expic = expic+getlitwinprob2by2(p0,p0e,nc,s,m,'gam')
          }
          expic = expic/np0e
          
          search_result = c(p0,p1,nc,ne,s,m,
                            round(alpL[iam], 5), round(alp1L[iam], 5), round(alp2L[iam], 5),
                            round(alphamaxvalue, 5), round(pwrL[iam], 5),round(betaL[iam], 5),
                            round(eta, 5),round(gam, 5),round(nogoL[iam], 5), round(expic, 5),
                            nc+ne,ci)
  
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

#### experiment ####
litwin_1s = function(FILEPREFIX, alphamax, alpha1max, alpha2max, alphamaxmax, betamax, pwrmin, ci=0.9,k=1) {
  FILENAME = paste(FILEPREFIX, '_confidence',round(ci*100), '.csv', 
                   sep='')
  print(FILENAME)
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
                   # etamax=etamax,
                   # gammax=gammax,
                   # expicmax=expicmax,
                   ci=ci)
  
  # write data
  sink(paste('results/',FILENAME, sep=''), append=F)
  cat('\n')
  cat(FILENAME)
  cat('\n')
  cat(paste0(c('p0', 'p1', 'nc','ne','s','m',
               'alpha', 'alpha1', 'alpha2', 'alphamax', 'power','beta','eta', 'gam', 'nogo','expic','N','ci'),collapse = ","))
  closeAllConnections()
  
  Litwin_1s_2by2_results = c()
  
  for (i in 1:nrow(dat)) {
    print(i)
    tic()
    s1 = getlitwinfunc1s(dat[i,],k=k,FILENAME=FILENAME)
    print(s1)
    toc()
    Litwin_1s_2by2_results = rbind(Litwin_1s_2by2_results, s1)
  }
  
}

fileprefix = 'Litwin_alpha20_1s_alphaonly_1to1.csv'
litwin_1s(FILEPREFIX=fileprefix,
            alphamax=0.20, alpha1max=0.99, alpha2max=0.99,alphamaxmax=0.99,
            betamax=0.20, pwrmin=0.01,
            ci=0.90,
            k=1)

fileprefix = 'Litwin_alpha10_1s_alphaonly_1to1.csv'
litwin_1s(FILEPREFIX=fileprefix,
          alphamax=0.10, alpha1max=0.99, alpha2max=0.99,alphamaxmax=0.99,
          betamax=0.10, pwrmin=0.01,
          ci=0.90,
          k=1)

fileprefix = 'Litwin_alpha20_1s_alphamax30_1to1.csv'
litwin_1s(FILEPREFIX=fileprefix,
          alphamax=0.20, alpha1max=0.99, alpha2max=0.99,alphamaxmax=0.20,
          betamax=0.20, pwrmin=0.01,
          ci=0.30,
          k=1)

fileprefix = 'Litwin_alpha10_1s_alphamax30_1to1.csv'
litwin_1s(FILEPREFIX=fileprefix,
          alphamax=0.10, alpha1max=0.99, alpha2max=0.99,alphamaxmax=0.10,
          betamax=0.10, pwrmin=0.01,
          ci=0.30,
          k=1)


# Litwin_1s_2by2_results = as.data.frame(Litwin_1s_2by2_results)
# colnames(Litwin_1s_2by2_results) = c('p0', 'p1', 'nc','ne','s','m',
#                                      'alpha', 'alpha1', 'alpha2', 'alphamax', 'power','beta','eta', 'gam', 'nogo','expic','N','ci')
# rownames(Litwin_1s_2by2_results) = 1:nrow(Litwin_1s_2by2_results)

