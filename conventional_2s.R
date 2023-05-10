##########################################################################
#                      Conventional Method Two-Stage                     #
##########################################################################
# -----------------------------------------------------------------------#
# code reference: Litwin et al. 2007
# -----------------------------------------------------------------------#

library(dplyr)
library(tictoc)

pbinres=function(p0,p1,nc,k=1) {
  # return probability matrix of possible combinations of responses in the experiemental arm and
  # the control arm under the null and alternative hypothesis
  # nc: size of the control arm
  # k: randomization ratio (E:C), default is 1
  ne=nc*k # 1:1 randomization
  ema=matrix( dbinom(0:ne,ne,p1), ne+1, 1 ) # +1 because of #res starts from 0 >> n+1 possibility
  emn=matrix( dbinom(0:ne,ne,p0), ne+1, 1 ) 
  cm=matrix( dbinom(0:nc,nc,p0), 1, nc+1 )
  ecma=ema%*%cm # matrix of all possible combination of responses under alternative hypothesis
  ecmn=emn%*%cm # matrix of all possible combination of responses under null hypothesis
  return( list(ecma,ecmn ) )      #  note if x=pbinres(p0,p1,ne) then ecma=x[[1]]  ecmn=x[[2]]
  
}

getstg1=function(es0min,es1max,p0,p1,LL,k=1) {
  # return stage 1 design of p0,p1,nc1,m1,a1,es0,es1
  # es0min: minimum early stopping probability under the null hypothesis
  # es1max: maximum early stopping probability under the alternative hypothesis
  # k: randomization ratio (E:C), default is 1
  
  #  nc1 loop
  sol=0
  nc1=2
  while(nc1<99) {
    if(sol>=1) break
    nc1=nc1+1 
    ne1=k*nc1
    pec=pbinres(p0,p1,nc1,k)
    peca=pec[[1]]
    pecn=pec[[2]]
    
    es0max=0
    es1min=1
    maxm=min( which(pbinom(1:ne1,ne1,p1) > es1max) ) 
    if (k==2) {
      a1L = 0:-10
      
    } else {
      a1L= 5:-10
    }
    
    for ( m1 in 1:maxm) { 
      #for ( m1 in 1:nc1) { 
      for (a1 in a1L){
        # calculate es1   
        es1=0
        ec1L=which(peca>LL)
        for (ec1 in ec1L) {
          c1=floor((ec1-.01)/(ne1+1))+1       # must subtract .01 , when ec3=k*(ne2+1) it does not work # get c1 from column number of the matrix
          e1=ec1-(c1-1)*(ne1+1)
          del=e1-1-k*(c1-1)
          # del <= a1    declare futility
          if( del <= a1 || (e1-1)<m1 ) {                      
            es1=es1+peca[e1,c1]
          } 
        }        # end ec1 loop  for es1
        
        # calculate es0
        ec1L=which(pecn>LL)
        es0=0
        for (ec1 in ec1L) {
          c1=floor((ec1-.01)/(ne1+1))+1       # must subtract .01 , when ec3=k*(ne2+1) it does not work
          e1=ec1-(c1-1)*(ne1+1)
          del=e1-1-k*(c1-1)
          # del <= a1    declare futility
          if( del <= a1 || (e1-1)<m1 ) {
            es0=es0+pecn[e1,c1]
          }
        }        # end ec1 loop  es0
        
        #    below uses #1 for max es0 criterion or #2 for min es1 criterion  to select (a1,m1)
        #    if( es1 <=es1max  && es1 <  es1min && es0 > es0min )          #2    needs brace
        if( es1 <= es1max && es0 >  es0max && es0 > es0min )  {  
          #1   needs brace
          es1min=es1                                    #2  not needed for #1
          es0max=es0   
          print('found')#1   not needed for #2
          ret=c(p0,p1,nc1,m1,a1,es0,es1)
          sol=sol+1 }  
        
      }        #  end a1
    }        # end m1 loop
  }        #  end nc1/ne1
  #        cat("p0,p1,k=" ,p0,p1,k,"\n")
  return(ret)
}        #  end fu


getconventionalstg2 = function(p0,p1,nc1,nc2,m1,a1,sL,m2L,call,L=1e-9,k=1) {
  # return a probability matrix of shape (nsL, nmL), where nsL and nmL are length of input 
  # parameters vectors under specified response rates under the null and alternative hypothesis
  # sL: vector of integers
  # mL: vector of integers
  # nc: size of the control arm
  # k: randomization ratio (E:C), default is 1
  nsL = length(sL)
  nm2L = length(m2L)
  psm2 = array(0, c(nsL,nm2L))
  ne1 = nc1*k
  ne2 = nc2*k
  p2n = array(0, c(ne1+1, nc1+1))
  
  # first cohort
  pe1c1 = pbinres(p0,p1,nc1,k)
  if (p1-p0<0.001) { # under H0
    pec1 = pe1c1[[2]]
  } else {
    pec1 = pe1c1[[1]]
  }
  
  ec1L = which(pec1 > L)
  for (ec1 in ec1L) {
    c1=floor((ec1-.01)/(ne1+1))+1       # must subtract .01 , when ec3=k*(ne2+1) it does not work
    e1=ec1-(c1-1)*(ne1+1)
    z1=e1-1-k*(c1-1)
    if(z1>a1 && (e1-1)>=m1) {
      p2n[e1,c1]=p2n[e1,c1]+pec1[e1,c1]
      }
  } # end of ec1
  
  # second cohort
  pe3c3 = pbinres(p0,p1,nc2,k)
  if ((p1-p0) < 0.001) { # under H0
    pec3 = pe3c3[[2]]
  } else {
    pec3 = pe3c3[[1]]
  }
  ec3L = which(pec3 > L) # reduce unnecessary computation
  ec2L = which(p2n > L) #proceed to stage 2
  for (ec3 in ec3L) {
    c3 = floor((ec3-0.01)/(ne2+1)) + 1
    e3 = ec3 - (c3-1)*(ne2+1)
    
    for (ec2 in ec2L) {
      c2=floor((ec2-.01)/(ne1+1))+1       # must subtract .01 , when ec3=k*(ne2+1) it does not work
      e2=ec2-(c2-1)*(ne1+1)               # (e2,c2) are first cohort that "continue"
      psav=pec3[e3,c3]*p2n[e2,c2]      # g function
      E = e2+e3-2 
      Z = e2+e3-2-k*(c2+c3-2) 
      if (call=='power' || call=='alpha') {
        psm2 = psm2 + outer(Z>=sL,E>=m2L) * psav
      }
      if (call=='type2' || call=='nogo') {
        psm2 = psm2 + outer(Z<=sL, 999>m2L) * psav 
      } else if (call == 'eta' || call=='gam') {
        psm2 = psm2 + outer(Z>=sL, E<m2L)*psav
      } 
    }
  }
  return(psm2) # 2 dimension s*m2
}

####---------------------Main Function----------------------####
getTDRfunc2s2by2 = function(x, k=1, FILENAME=FILENAME) {
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
  es0min = x['es0min'][[1]]
  es1max = x['es1max'][[1]]
  ci = x['ci'][[1]]
  
  p01=p0+.05     
  p11=p1+.05      ;  if(p11>.999) p11=1
  p02=(p0+p1)/2
  p12=p1+(p1-p0)/2  ; if(p12>.999) p12=1
  
  L=1e-9
  LL=1e-7
  
  #  get nc1, a1, m1 
  s1=getstg1(es0min,es1max,p0,p1,LL,k) 
  print(s1)
  nc1=s1[3]
  ne1=k*nc1
  m1=s1[4]
  a1=s1[5]
  es0=s1[6]
  es1=s1[7]
  
  #first cohort:   get p2n & p2a
  p2n=matrix(0,nrow=ne1+1,ncol=ne1+1)  #p2n[e2,c2]=Pr{e2,c2|take 2nd cohart H0}
  p2a=matrix(0,nrow=ne1+1,ncol=ne1+1)  #                                    H1  
  
  pe1vn=dbinom(0:ne1,ne1,p0)
  pe1va=dbinom(0:ne1,ne1,p1)
  pc1vn=dbinom(0:nc1,nc1,p0)
  
  
  for (e1 in 0:ne1 ) {
    for (c1 in 0:nc1 ) {  
      z1=e1-k*c1
      if(z1 > a1 && e1>=m1 )   # proceed to stage 2               
      {
        p2n[e1+1,c1+1]=p2n[e1+1,c1+1] + pe1vn[e1+1]*pc1vn[c1+1]             
        p2a[e1+1,c1+1]=p2a[e1+1,c1+1] + pe1va[e1+1]*pc1vn[c1+1] 
      } 
    }  # end c1
  }  # end e1
  
  # second cohort
  soln = 0
  nc2 = round(nc1*0.2) # 4
  nc2s = 0
  nc2max = 100
  
  while (nc2 <= nc2max) { # nc2 loop step by 4
    
    if( soln==2 ) { nc2 = nc2+1   }     # first soluton found and no solution found for ne2 
    
    ne2 = k*nc2 
    print(ne2)
    
    NE = ne1 + ne2
    NC = nc1 + nc2
    se0 = sqrt(p0*(1-p0)*NC)
    se1 = sqrt(p1*(1-p1)*NE)
    se10 = sqrt(p1*(1-p1)*NE+p0*(1-p0)*NC)
    
    smin = max(-NC, floor(p0*(NE-k*NC)))
    smax = min(NE, ceiling(p1*NE-p0*k*NC+4*se10))
    
    m2min = floor(p0*NE)
    m2max = min(NE, ceiling(NE*p1+4*se1))
    
    sL = smin:smax
    m2L = m2min:m2max
    
    p0L = max(0.001, p0-qnorm(0.5+0.5*ci)*sqrt(p0*(1-p0)/((1+k)*(NC))))
    p0U = min(0.999, p0+qnorm(0.5+0.5*ci)*sqrt(p0*(1-p0)/((1+k)*(NC))))
    p0.ci = seq(p0L, p0U, length.out=11)
    alpha.actual = rep(0,11) # search worst case scenario for type I error
    pwr.actual = rep(0,11)
    
    # control alpha
    alpL = getconventionalstg2(p0,p0,nc1,nc2,m1,a1,sL,m2L,'alpha',LL,k)
    alpLidx = which(alpL<=alphamax)
    
    alpLsi = which(alpL <= alphamax, arr.ind = T)
    constrains = sL[alpLsi[,1]]
    constrainm2 = sL[alpLsi[,2]]
    constrain = constrains < constrainm2
    alpLsicontrain = alpLsi[constrain,]
    nm2 = length(m2L)
    ns = length(sL)
    if (length(alpLsicontrain)==2) {
      constrainx = (alpLsicontrain[2]-1)*ns + alpLsicontrain[1]
    } else {
      constrainx = (alpLsicontrain[,2]-1)*ns + alpLsicontrain[,1]
    }
    xd = intersect(alpLidx, constrainx)
    
    if (length(xd) > 0) {  # control beta
      betaL = getconventionalstg2(p0,p1,nc1,nc2,m1,a1,sL,m2L,'type2',LL,k) 
      betaLidx = which(betaL<=betamax)
      xd=intersect(xd, betaLidx)
      
      if (length(xd) > 0) { # control alpha1
        alp1L = getconventionalstg2(p01,p01,nc1,nc2,m1,a1,sL,m2L,'alpha',LL,k)
        alp1Lidx = which(alp1L<=alpha1max)
        xd = intersect(xd, alp1Lidx)
        
        if (length(xd) > 0) { # control alpha2
          alp2L = getconventionalstg2(p02,p02,nc1,nc2,m1,a1,sL,m2L,'alpha',LL,k)
          alp2Lidx = which(alp2L<=alpha2max)
          xd = intersect(xd, alp2Lidx)
          
          if (length(xd) > 0) { # control alpha max
            for (p0s in p0.ci) { 
              alpmaxL = getconventionalstg2(p0s,p0s,nc1,nc2,m1,a1,sL,m2L,'alpha',LL,k)
              alpmaxLidx = which(alpmaxL<=alphamaxmax)
              xd=intersect(xd,alpmaxLidx)
            }
            if (length(xd) > 0) { # control pwr
              pwrL = getconventionalstg2(p0,p1,nc1,nc2,m1,a1,sL,m2L,'power',LL,k) # probability list
              pwrLidx = which(pwrL >= pwrmin) # index
              xd = intersect(pwrLidx, xd)
              
              if (length(xd) > 0) { # control eta
                etaL = getconventionalstg2(p0,p0,nc1,nc2,m1,a1,sL,m2L,'eta',LL,k)
                etaLidx = which(etaL<=etamax)
                xd=intersect(xd, etaLidx)
                
                if (length(xd) > 0) {
                  gamL = getconventionalstg2(p0,p1,nc1,nc2,m1,a1,sL,m2L,'gam',LL,k) 
                  gamLidx = which(gamL<=gammax)
                  xd=intersect(xd, gamLidx)
                  
                  if (length(xd) > 0) {
                    
                    expicL = (etaL + gamL)/2
                    expicLidx = which(expicL<=expicmax)
                    xd=intersect(expicLidx,xd)
                    
                    nogoL = getconventionalstg2(p0,p0,nc1,nc2,m1,a1,sL,m2L,'nogo',LL,k)
                    nogoLidx = which(nogoL>=nogomin)
                    xd=intersect(xd, nogoLidx)
                    # print(c(max(pwrL[xd]), min(alpL[xd]), min(betaL[xd]),max(nogoL[xd])))
                    
                    if (length(xd) > 0) soln = soln + 1
                  } # expic nogo
                } # gam
              } # eta
            } # pwr
          } # control alphamax --> variations of p0 influence the type I error mostly
        } # control alpha2
      } # control alpha1
    } # control beta
    En = ne1+nc1 + (1-es0)*(ne2+nc2)
    
    if( soln==0 ) { nc2=nc2+8 }  
    
    if (soln == 1) {
      nc2s = nc2
      nc2 = max(1,nc2-9)
      soln = 2
    }
    
    if (soln==3 || nc2==nc2s) {
      jz = which(alpL==min(alpL[xd])) 
      jz = intersect(xd, jz)
      iam = jz[1]
      if (length(iam) <= 0) {
        iam = xd[1]
      }
      if (length(iam) > 0) {
        m2idx = floor((iam-0.01)/ns) + 1
        sidx = iam - (m2idx-1)*ns
        s = sL[sidx]
        m2 = m2L[m2idx]
        
        N1 = ne1+nc1
        N2 = N1 + ne2 + nc2
        print('minmax solution found')
        
        alphamaxvalue = -0.01
        for (p0s in p0.ci) { 
          alpmaxL = getconventionalstg2(p0s,p0s,nc1,nc2,m1,a1,s,m2,'alpha',LL,k)
          alphamaxvalue = max(alphamaxvalue, alpmaxL[1])
        }
        
        search_result = c(p0,p1,nc1,ne1,a1,m1,nc2,ne2,s,m2,
                          round(alpL[iam], 5), round(alp1L[iam], 5), round(alp2L[iam], 5),
                          round(alphamaxvalue, 5), round(pwrL[iam], 5),round(betaL[iam], 5),
                          round(etaL[iam], 5),round(gamL[iam], 5), 
                          round(nogoL[iam], 5), round(expicL[iam], 5), round(es0, 5), round(es1,5),
                          En, N1, N2, etamax, gammax, expicmax,ci)
        print(search_result)
        
        sink(paste('results/',FILENAME, sep=''), append=T)
        cat('\n')
        cat(paste0(search_result, collapse = ","))
        closeAllConnections()
        
        return(search_result)
        break
    }
    } # minmax soln found
  } # while loop
  
}

####---------------------Experiment---------------------####
tdr_2s_2by2 = function(eta,gam,expic, ci=0.3) {
  FILENAME = paste('TDR_2s_2by2_alphamax_1to1_','_eta', round(eta*100),'_gam', round(gam*100), 
                   '_ic',round(expic*100), '_confidence',round(ci*100), '.csv', 
                   sep='')
  
  # TODO: add eta1,2,gam1,2 control above
  
  tdr_2s_2by2_results = c()
  print(FILENAME)
  K=1
  alpha1max = 0.99
  alpha2max = 0.99
  alphamaxmax = 0.20
  pwrmin = 0.01
  nogomin =0.01
  etamax = eta
  gammax = gam
  expicmax = expic
  ci = ci
  dat1 = data.frame(p0=c(5,5,5,10,10,15,15,20,20,25,25,30,30)/100, # 13 cases
                    p1=c(15,20,25,25,30,30,35,35,40,40,45,45,50)/100,
                    es0min=0.50,
                    es1max=0.05,
                    alphamax=0.20,
                    betamax=0.20,
                    alpha1max=alpha1max,
                    alpha2max=alpha2max,
                    alphamaxmax=alphamaxmax,
                    pwrmin=pwrmin,
                    nogomin=nogomin,
                    etamax=etamax,
                    gammax=gammax,
                    expicmax=expicmax,
                    ci=ci)
  
  dat2 = data.frame(p0=c(rep(seq(0.35,0.55,0.05), each=2)), # 10 cases
                    p1=c(0.50,rep(seq(0.55,0.70,0.05), each=2),0.75),
                    es0min=0.50,
                    es1max=0.05,
                    alphamax=0.20,
                    betamax=0.20,
                    alpha1max=alpha1max,
                    alpha2max=alpha2max,
                    alphamaxmax=alphamaxmax,
                    pwrmin=pwrmin,
                    nogomin=nogomin,
                    etamax=etamax,
                    gammax=gammax,
                    expicmax=expicmax,
                    ci=ci)
  
  dat3 = data.frame(p0=c(rep(seq(0.60,0.75,0.05), each=2),0.80, 0.85), # 10 cases
                    p1=c(0.75,rep(seq(0.80,0.95,0.05),each=2), 0.95),
                    es0min=0.50,
                    es1max=0.05,
                    alphamax=0.20,
                    betamax=0.20,
                    alpha1max=alpha1max,
                    alpha2max=alpha2max,
                    alphamaxmax=alphamaxmax,
                    pwrmin=pwrmin,
                    nogomin=nogomin,
                    etamax=etamax,
                    gammax=gammax,
                    expicmax=expicmax,
                    ci=ci)
  
  dat = rbind(dat1,dat2, dat3)
  
  #### write data ####
  sink(paste('results/',FILENAME, sep=''), append=F)
  cat('\n')
  cat(FILENAME)
  cat('\n')
  cat(paste0(c('p0', 'p1', 'nc1','ne1','a1','m1','nc2','ne2','s','m2',
               'alpha', 'alpha1', 'alpha2', 'alphamax','power','beta',
               'eta', 'gam', 'nogo','expic','es0','es1',
               'En','N1','N2','etamax','gammax','expicmax', 'ci'),collapse = ","))
  closeAllConnections()
  
  TDR_2s_results = c()
  
  for (i in 1:nrow(dat)) {
    print(i)
    tic()
    s2 = getTDRfunc2s2by2(dat[i,],k=K, FILENAME=FILENAME)
    toc()
    # zTDR_2s_results = rbind(TDR_2s_results, s2)
  }
  
}

####-------------------------Parameter Search------------------------####
for (eta in c(0.99,0.1,0.2,0.3,0.4)) {
  tdr_2s_2by2(eta,0.99,0.99)
}

for (gam in c(0.99,seq(0.01,0.2,0.01),0.3)) {
  tdr_2s_2by2(0.99,gam,0.99)
}

for (expic in c(0.99,0.05,0.1,0.2,0.3,0.4)) {
  tdr_2s_2by2(0.99,0.99,expic)
}

for (expic in c(0.1,0.2,0.3)) {
  for (gam in c(seq(0.02,0.2,0.02))) {
    tdr_2s_2by2(0.99,gam,expic)
  }
}

##########################################################################
#                                  End                                   #
##########################################################################
