library(dplyr)
library(tictoc)

pbinres=function(p0,p1,nc,k=1) {
  ne=nc*k # 1:1 randomization
  ecmn=matrix(0,ne+1,nc+1)
  ecma=matrix(0,ne+1,nc+1)
  ema=matrix( dbinom(0:ne,ne,p1), ne+1, 1 ) # +1 because of #res starts from 0 >> n+1 possibility
  emn=matrix( dbinom(0:ne,ne,p0), ne+1, 1 ) 
  cm=matrix( dbinom(0:nc,nc,p0), 1, nc+1 )
  ecma=ema%*%cm # matrix of all possible combination of responses under alternative hypothesis
  ecmn=emn%*%cm # matrix of all possible combination of responses under null hypothesis
  return( list(ecma,ecmn ) )      #  note if x=pbinres(p0,p1,ne) then ecma=x[[1]]  ecmn=x[[2]]
  
}

# 1:1 randomization
# litwin rpgmnum1.R
getlitwinstg1prob=function(es0min,es1max,p0,p1,LL,k=1) {
  #  nc1 loop
  sol=0
  nc1=2
  while(nc1<99) {
    print(nc1)
    if(sol>=1) break
    nc1=nc1+1  ;  ne1=k*nc1  ; ne11=ne1+1  
    pec=pbinres(p0,p1,nc1,k)
    peca=pec[[1]]
    pecn=pec[[2]]
    
    es0max=0
    es1min=1
    maxm=min(which(pbinom(1:ne1,ne1,p1) > es1max) ) 
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
          c1=floor((ec1-.01)/ne11)+1       # must subtract .01 , when ec3=k*ne21 it does not work # get c1 from column number of the matrix
          e1=ec1-(c1-1)*ne11
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
          c1=floor((ec1-.01)/ne11)+1       # must subtract .01 , when ec3=k*ne21 it does not work
          e1=ec1-(c1-1)*ne11
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

# getlitwinstg1prob(0.50,0.05,0.05,0.15,1e-7,2)
# getlitwinstg2prob(0.10,0.10,11,6,2,-5,0:14,1:6,'alpha',k=2)
# get stage 2 metrics
getlitwinstg2prob = function(p0,p1,nc1,nc2,m1,a1,b2L,m2L,call,L=1e-9,k=1) {
  
  nb2L = length(b2L)
  nm2L = length(m2L)
  pm2b2 = array(0, c(nm2L,nb2L))
  ne1 = nc1*k
  ne2 = nc2*k
  p2n = array(0, c(ne1+1, nc1+1))
  
  # first cohort
  pe1c1 = pbinres(p0,p1,nc1,k)
  if ((p1-p0)<.001) { # under H0
    pec1 = pe1c1[[2]]
  } else {
    pec1 = pe1c1[[1]]
  }
  
  ec1L = which(pec1 > L)
  ne11 = ne1 + 1
  for (ec1 in ec1L) {
    c1=floor((ec1-.01)/ne11)+1       # must subtract .01 , when ec3=k*ne21 it does not work
    e1=ec1-(c1-1)*ne11
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
  ne21 = ne2 + 1
  ec2L = which(p2n > L) #proceed to stage 2
  for (ec3 in ec3L) {
    c3 = floor((ec3-0.01)/ne21) + 1
    e3 = ec3 - (c3-1)*ne21
    for (ec2 in ec2L) {
      c2=floor((ec2-.01)/ne11)+1       # must subtract .01 , when ec3=k*ne21 it does not work
      e2=ec2-(c2-1)*ne11               # (e2,c2) are first cohort that "continue"
      psav=pec3[e3,c3]*p2n[e2,c2]      # g function
      E = e2+e3-2 
      Z = E-k*(c2+c3-2)
      if (call=='alpha'|| call == 'power') {
        pm2b2 = pm2b2 + outer(E>=m2L,Z>=b2L) * psav
      }
      if (call=='type2' || call=='nogo') {
        Imat = array(1,c(length(m2L), length(b2L)))
        pm2b2 = pm2b2 + (Imat - outer(E>=m2L,Z>=b2L)) * psav
      }
      
      # a 3-dim binary matrix multiply with g(px,px)
    }
  }
  return(pm2b2) # 3 dimension b1*b2*m2
}

getlitwinfunc2s = function(x, k=1, FILENAME=FILENAME) {
  p0=x['p0'][[1]]
  p1=x['p1'][[1]]
  es0min = x['es0min'][[1]]
  es1max = x['es1max'][[1]]
  alphamax = x['alphamax'][[1]]
  pwrmin = x['pwrmin'][[1]]
  alpha1max = x['alpha1max'][[1]]
  alpha2max = x['alpha2max'][[1]]
  alphamaxmax = x['alphamaxmax'][[1]]
  betamax = x['betamax'][[1]]
  nogomin = x['nogomin'][[1]]
  ci = x['ci'][[1]]
  ifL=1               #  index for first loop
  # neL=1               #  never ending loop
  
  p01=p0+.05     
  p11=p1+.05      ;  if(p11>.999) p11=1
  p02=(p0+p1)/2
  p12=p1+(p1-p0)/2  ; if(p12>.999) p12=1
  
  L=1e-9
  LL=1e-7 
  
  #  get nc1, a1, m1 
  s1=getlitwinstg1prob(es0min,es1max,p0,p1,LL,k) 
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
  soln=0                  #  if soln=0  no solution found
  nc2s=0           
  nc2=4 
  nc2max=130
  while (nc2 <= nc2max ) { # nc2 loop step by 4
    
    #cat("ne2,soln=",ne2,soln,"\n")
    if( soln==2 ) { nc2 = nc2+1   }     # first soluton found and no solution found for ne2 
    # if(soln==3 || nc2==nc2s ) break
    
    ne2=nc2*k
    print(ne2)
    
    NC = nc1 + nc2
    NE = ne1 + ne2
    se10 = sqrt(p1*(1-p1)*NE+p0*(1-p0)*NC)
    
    m2min = floor(p0*NE)
    m2max = min(ne1+ne2, ceiling(NE*p0+4*sqrt(NE*p0*(1-p0))))
    
    b2min = max(-NC, floor(p0*(NE-k*NC)))
    b2max = min(NE, ceiling(p1*NE-p0*k*NC + 4*se10))
    
    b2L = 1:12# b2min:b2max
    m2L = m2min:m2max
    
    
    p0L = max(0.001, p0-qnorm(0.5+0.5*ci)*sqrt(p0*(1-p0)/((1+k)*(NC))))
    p0U = min(0.999, p0+qnorm(0.5+0.5*ci)*sqrt(p0*(1-p0)/((1+k)*(NC))))
    p0.ci = seq(p0L, p0U, length.out=11)
    alpha.actual = rep(0,11) # search worst case scenario for type I error
    
    # control power
    pwrL = getlitwinstg2prob(p0,p1,nc1,nc2,m1,a1,b2L,m2L, 'power',L,k) # probability list
    pwrLidx = which(pwrL >= pwrmin) # index
    
    # pwrLsi = which(pwrL >= pwrmin, arr.ind = T)
    # constrainm2 = m2L[pwrLsi[,1]]
    # constrainb2 = b2L[pwrLsi[,2]]
    # 
    # constrain = constrainb2 < constrainm2
    # pwrLsicontrain = pwrLsi[constrain,]
    # nm2 = length(m2L)
    # nb2 = length(b2L)
    # if (length(pwrLsicontrain)==2) {
    #   constrainx = (pwrLsicontrain[2]-1)*nm2 + pwrLsicontrain[1]
    # } else {
    #   constrainx = (pwrLsicontrain[,2]-1)*nm2 + pwrLsicontrain[,1]
    # }
    # xd = intersect(pwrLidx, constrainx)
    xd = pwrLidx
    # print(max(pwrL))
    if (length(pwrLidx) > 0) {  # control alpha
      alpL = getlitwinstg2prob(p0,p0,nc1,nc2,m1,a1,b2L,m2L,'alpha',L,k)
      alpLidx = which(alpL<=alphamax)
      xd=intersect(xd, alpLidx)
      
      if (length(xd) > 0) { # control alpha1
        alp1L = getlitwinstg2prob(p01,p01,nc1,nc2,m1,a1,b2L,m2L,'alpha',L,k)
        alp1Lidx = which(alp1L <= alpha1max)
        xd = intersect(xd, alp1Lidx)
        
        if (length(xd)>0) {
          alp2L = getlitwinstg2prob(p02,p02,nc1,nc2,m1,a1,b2L,m2L,'alpha',L,k)
          alp2Lidx = which(alp2L <= alpha2max)
          xd = intersect(xd, alp2Lidx)
          
          if (length(xd > 0)) {
            for (p0s in p0.ci) { # control alpha max
              alpmaxL = getlitwinstg2prob(p0s,p0s,nc1,nc2,m1,a1,b2L,m2L,'alpha',L,k)
              alpmaxLidx = which(alpmaxL<=alphamaxmax)
              xd=intersect(alpmaxLidx,xd)
            }
            if (length(xd) > 0) {
              soln = soln+1
              betaL = getlitwinstg2prob(p0,p1,nc1,nc2,m1,a1,b2L,m2L,'type2',L,k)
              nogoL = getlitwinstg2prob(p0,p0,nc1,nc2,m1,a1,b2L,m2L,'nogo',L,k)
            }
          }
        }
      }
    } # control alpha
    En = (ne1+nc1) + (1-es0)*(ne2+nc2)
    
    if( soln==0 ) { nc2=nc2+8 }  
    
    if (soln == 1) {
      pwrLsi = which(pwrL >= pwrmin, arr.ind = T)
      jz=which(alpL==min(alpL[xd])) # within combinationa of 
      # b1,b2,m2 meeting criteria, select the smallest alpha
      iam=which(pwrLidx==jz)
      m2=m2L[pwrLsi[iam,1]]
      b2=b2L[pwrLsi[iam,2]]
      # neL == 0
      
      print(c(p0,p1,nc1,ne1,a1,m1,nc2,ne2,b2,m2,
              round(alpL[jz],5),round(alp1L[jz],5), round(alp2L[jz],5),round(pwrL[jz],5),round(betaL[jz],5),round(nogoL[jz],5),
              En))
      nc2s=nc2
      nc2=max(1,nc2-9)
      soln=2
      }
    if (soln==3 || nc2==nc2s) {
      
      pwrLsi = which(pwrL >= pwrmin, arr.ind = T)
      jz=which(alpL==min(alpL[xd])) # within combinationa of 
      # b1,b2,m2 meeting criteria, select the smallest alpha
      jz=jz[1]
      iam=which(pwrLidx==jz)
      m2=m2L[pwrLsi[iam,1]]
      b2=b2L[pwrLsi[iam,2]]
      
      print('minmax solution found')
      search_result = c(p0,p1,nc1,ne1,a1,m1,nc2,ne2,b2,m2,
                        round(alpL[jz],5),round(alp1L[jz],5), round(alp2L[jz],5),round(pwrL[jz],5),round(betaL[jz],5),round(nogoL[jz],5),
                        En)
      print(search_result)
      return(search_result)
      break
        }
    
  }
}

#### experiment ####
litwin_2s = function(FILEPREFIX, alphamax, alpha1max, alpha2max, alphamaxmax, betamax, pwrmin, ci=0.9,k=1) {
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
  litwin_results = c()
  for (i in 1:nrow(dat)) {
    print(i)
    tic()
    s2 = getlitwinfunc2s(dat[i,],k=k)
    toc()
    litwin_results = rbind(litwin_results, s2)
  }
  
  litwin_results = as.data.frame(litwin_results)
  colnames(litwin_results) = c('p0','p1','nc1','ne1','a1','m1','nc2','ne2','b2','m2',
                               'alpha', 'alpha1','alpha2', 'power','beta','nogo',
                               'En')
  rownames(litwin_results) = 1:nrow(litwin_results)
  litwin_results$N1 = litwin_results$nc1+litwin_results$ne1
  litwin_results$N2 = litwin_results$nc2+litwin_results$ne2 + litwin_results$N1
  litwin_results$ci = ci
  write.csv(litwin_results, paste('results/',FILENAME,sep=''))
  
}

fileprefix = 'Litwin_alpha20_2s_alphaonly_1to1.csv'
litwin_2s(FILEPREFIX=fileprefix, 
          alphamax=0.20, alpha1max=0.99, alpha2max=0.99,alphamaxmax=0.99, 
          betamax=0.20, pwrmin=0.80,
          ci=0.90,
          k=1)

fileprefix = 'Litwin_alpha10_2s_alphaonly_1to1.csv'
litwin_2s(FILEPREFIX=fileprefix, 
          alphamax=0.10, alpha1max=0.99, alpha2max=0.99,alphamaxmax=0.99, 
          betamax=0.10, pwrmin=0.90,
          ci=0.90,
          k=1)

fileprefix = 'Litwin_alpha20_2s_alphamax30_1to1.csv'
litwin_2s(FILEPREFIX=fileprefix, 
          alphamax=0.20, alpha1max=0.99, alpha2max=0.99,alphamaxmax=0.20, 
          betamax=0.20, pwrmin=0.80,
          ci=0.30,
          k=1)

fileprefix = 'Litwin_alpha10_2s_alphamax30_1to1.csv'
litwin_2s(FILEPREFIX=fileprefix, 
          alphamax=0.10, alpha1max=0.99, alpha2max=0.99,alphamaxmax=0.10, 
          betamax=0.10, pwrmin=0.90,
          ci=0.30,
          k=1)




