library(meanShiftR)

## This code performs all the analysis to obtain the #genotypes, #clusters, Tajima's D, #fixations, genetic and phenotypic diversity at the end of the simulations.
## It outputs one data matrix per condition. 


### input file name ###
InFile="Results_Data.RData"

#The data come as iteration of 1 parameter: mutation rate
load(InFile)

### Output file name ###
OutFile="Summary_Data.csv"

i=0
parameter_comb=list()
for(u in Us){
  i=i+1
  parameter_comb[[i]]=c(R,u,pleiotropy)
}

save(allAv_alphas, file = "Alpha_Dynamics.RData") ##this file will contain the data for the dynamics of the population trait sum as in Fig.2E
save(allStrains, file = "M_Dynamics.RData") ##this file will contain the data for the M dynamics as in Fig.3A

#### Summary data as a matrix with (sigma,pleiotropy,U,N,M*,Clusters,Pi(G),Pi(P),Tajima's D and Fixations) ####
tot_simulations=simulations*length(Us)
Data=list()
for(u in 1:length(Us)){
  Data[[u]]=data.frame(rep(Sigma,simulations),rep(pleiotropy,simulations),rep(Us[u],simulations),rep(Pop,simulations),rep(0,simulations),rep(0,simulations),rep(0,simulations),rep(0,simulations),rep(0,simulations),rep(0,simulations))
  colnames(Data[[u]])=c("sigma","pleiotropy","U","N","M*","Clusters","Pi(G)","Pi(P)","Tajima's D","Fixations")
}

### Number of genotypes ###
column=which(colnames(Data[[1]])=="M*")
time=t_tot
time_i=9000 ### starting point to compute M*
time_f=t_tot ## final point to compute M*
for (i in 1:counter){
  R=parameter_comb[[i]][1]
  u=parameter_comb[[i]][2]
  
  Strains=allStrains[[i]]
  #History=allHistory[[i]]
  
  Mstar=vector()
  for (replica in 1:simulations){
    Mstar=c(Mstar,mean(Strains[replica,time_i:time_f]))
  }
  Data[[i]][,column]=Mstar
}


### Number of clusters ###
column=which(colnames(Data[[1]])=="Clusters")
time=t_tot
for (i in 1:counter){
  R=parameter_comb[[i]][1]
  u=parameter_comb[[i]][2]
  pleiotropy=parameter_comb[[i]][3]
  Delta=Sigma*sqrt(1+(R-1)*pleiotropy^2)
  h=rep(Sigma,R)
  epsiloncluster=Delta
  
  Alphas=allAlphas[[i]]
  Surv=allSurv[[i]]
  #History=allHistory[[i]]
  
  Clusters=vector()
  for (replica in 1:simulations){
    surv=Surv[[replica]]
    alphas=Alphas[[replica]][surv,]
    m=matrix(rep(0,R),ncol=R,nrow = 1)
    data=as.matrix(alphas)
    cluster=meanShift(data,data,bandwidth = h,epsilonCluster = epsiloncluster)
    Clusters=c(Clusters,max(cluster$assignment))
  }
  if(epsiloncluster==0){Clusters=rep(1,length(Clusters))}
  Data[[i]][,column]=Clusters
  
}

### Tajima's D ###
columnT=which(colnames(Data[[1]])=="Tajima's D")
columnF=which(colnames(Data[[1]])=="Fixations")
n=100 ## sampled population of size 100

for(u in 1:length(Us)){
  a1=0
  a2=0
  for(i in 1:(n-1)){
    ### Compute the parameters as in Tajima 1989 ###
    a1=a1+1/i
    a2=a2+1/(i^2)}
    b1=(n+1)/(3*(n-1))
    b2=2*(n^2+n+3)/(9*n*(n-1))
    e1=(b1-(1/a1))/a1
    e2=(b2-(n+2)/(a1*n)+(a2/(a1^2)))/((a1^2)+a2)
    
    TajimasD=vector()
    fixations=vector()
  for( ex in 1:simulations){
    History=allHistory[[u]][[ex]]
    densities=allDyn[[u]][[ex]][[t_tot]][2,]
    
    surv=allSurv[[u]][[ex]]
    set.seed(ex)
    types=sample(rep(surv,round(densities)),n) ##sampled population
    count=0
    Pi_dist=0
    Union=unlist(History[[types[n]]])
    fix=1:max(unlist(History))
    
    ## pairwise genetic distance
    for (i in 1:(length(types)-1)){
      a=unlist(History[[types[i]]])
      for(j in (i+1):length(types)){
        count=count+1
        
        b=unlist(History[[types[j]]])
        fix=intersect(fix,intersect(a,b))
        Pi_dist=Pi_dist+length(setdiff(union(a,b),intersect(a,b)))
      }
      Union=union(Union,a)
      
    }
    Pi_dist=Pi_dist/choose(n,2)
    fixations=c(fixations,length(fix)-1)
    S=length(setdiff(Union,fix))
    S_dist=S/a1
    if(S!=0){TajimasD=c(TajimasD,(Pi_dist-S_dist)/(sqrt(e1*S+e2*S*(S-1))))}else{TajimasD=c(TajimasD,NA)}
    
  }
  Data[[u]][,columnT]=TajimasD
  Data[[u]][,columnF]=fixations
}


### Genetic and Phenotypic diversity ###
n=100 ## size of the sampled population
columnG=which(colnames(Data[[1]])=="Pi(G)")
columnP=which(colnames(Data[[1]])=="Pi(P)")
for(u in 1:length(Us)){
  GD=c()
  PD=c()
  for( ex in 1:simulations){
    History=allHistory[[u]][[ex]]
    densities=allDyn[[u]][[ex]][[t_tot]][2,]
    
    surv=allSurv[[u]][[ex]]
    
    set.seed(ex) ##To have the same sampling as for the Tajima
    types=sample(rep(surv,round(densities)),n)
    alphas=allAlphas[[u]][[ex]][types,]
    
    Pi_dist=0
    
    fix=1:max(unlist(History))
    for (i in 1:(length(types)-1)){
      for(j in (i+1):length(types)){
        a=unlist(History[[types[i]]])
        b=unlist(History[[types[j]]])
        fix=intersect(fix,intersect(a,b))
        pi_dist=length(setdiff(union(a,b),intersect(a,b)))
        Pi_dist=Pi_dist+pi_dist
        
      }
      
    }
    GD=c(GD,Pi_dist/choose(n,2))
    PD=c(PD,sum(dist(alphas))/choose(n,2))
    
  }
  Data[[u]][,columnG]=GD
  Data[[u]][,columnP]=PD/Sigma
}

for(u in 1:length(Us)){
  data=Data[[u]]
  write.csv(data,paste("U",Us[u],OutFile,sep=""))
}
