## HIGH NU ##
#!/usr/bin/env Rscript
#options(repos=structure(c(CRAN="http://cran.mirror.garr.it/mirrors/CRAN/")))

### This algorithm models adaptation across different mutation rates (Us) and keeps the other parameters constant

######## Algorithm Parameters #######

simulations=2 ## number of simulations
t_tot=10000 ##Total time/generations
R=2 ##number of resources (>=2)
d=1 ##death rate
Pop=10^7    # Population size

### Ancestor phenotypes
inFit=0.1 ## set the sum of the alphas

#### Evolution Parameters ####
## mutations ##
Us=c(10^(-8),10^(-7))   ## Mutation Rate vector, the algorithm iterates the different values and save in a list
Sigma=0.05   ## Mutation step standard deviation
pleiotropy=0 ## Pleiotropic effect
n0=1 ## initial density of the mutants
E_constraint=1 ##Energy constraint
#########

### output file name ###
OutFile="Results_Data.RData"

#########

print(paste("#simulations:",simulations,", Population size:",Pop,", Mutation rate (U):",Us,", Mutation step:",Sigma,", Pleiotropy: ",pleiotropy))
print(paste("Number of resources:",R))


##### Start the adaptation simulations ####

##Initialize Data structures
My_data=c("Av_alphas","Strains") ## These dataframes will store the average phenotype and the number of strains


allCPtime=list() ## save computation times at t_tot
allStrains=list() ## save number of strains over time
allAlphas=list() ## save the alphas 
allDyn=list() ## save all the densities over time (dynamics)
allHistory=list() ## save all the philogenetic info
allSurv=list() ## save surviving strains at t_tot
allAv_alphas=list() ## save mean sum(alphas) over time

counter=0 ## parameter iteration

## iterations over different resource numbers
############ resource supply ##########
sS=unlist(lapply(1:R, function(x)(Pop/R))) ##Here, resources are equally distributed
sS_norm=sS/(sum(sS))

### Computational time table
CPtime_table=data.frame(matrix(rep(0,length(Us)),nrow =1,ncol = length(Us)))
rownames(CPtime_table)=pleiotropy
colnames(CPtime_table)=Us
  
print(paste("##########################################",R," Resources ###########################################"))
  
### iterations over different mutation rates
for (u in 1:length(Us)){
    U=Us[u]
    pt0=proc.time()[3] ## reset the CP time
    counter=counter+1 ## update the count of parameter combinations
      
    print(paste("With U:",Us[u],", sigma ",Sigma," and pleiotropy:",pleiotropy,""))
    
    ####################### Different simulations  ####################
    Alphas=list()
    Dyn=list()
    History=list()
    Surv=list()
    
    ## iterate independent populations
    for(replica in 1:simulations){
      print(paste("simulation",replica, " has started!"))
      Dout=list()
      
      #### Initial conditions 
      N=1  ### number of initial species
      Ndym=c(N)  ## strain diversity over time
      
      #### alphas (phenotypes) ####
      # build an empty data frame for the phenotypes
      alphas=data.frame(matrix(rep(0,N*R),nrow=N,ncol=R))
      #assign th ancestor phenotype
      ancestor=rep(inFit/R,R)
      alphas[1,]=ancestor
        
      history=list()
      history[[1]]=c(1) ## store the phylogenetic info
        
      Ns=c(1) ## count all the genotypes
      trackID=c(1) ## to keep track of the genotypes' identity
        
      State=c(Pop,rep(0,N-1)) ## Initial density
      enzymes=unlist(lapply(1:R,function(x)(ancestor[x]*State[1])))
        
      Dout[[1]]=rbind(trackID,State) ## matrix of the dynamics
      av_alphas=c(sum(ancestor)) ## sum(alphas) over time
      
      ### Iterate over time/generations
      for (generation in 2:t_tot){
        t_alphas=alphas[trackID,] ## get the alphas only of the alive strains
        Nt=Ns 
        parents=vector()
        ## mutation iteration for each genotype
        for(x in 1:N){
          parental=alphas[trackID[x],] ## parental phenotype
          ###### Poisson process ######
          exp_mut=rpois(1,State[x]*U) ## number of mutations on genotype x
          
          ## iteration of mutations on genotype x  
          if(exp_mut>0){
            parents=c(parents,rep(trackID[x],exp_mut))
            for (trial in 1:exp_mut){
              mutant=rep(100,R)
              ## draw the mutation effect within the allowed region (sum(alphas)<=E)
              while(sum(mutant)>E_constraint){
                deltas=rnorm(R,mean=0,sd=Sigma*pleiotropy)
                #deltas=rep(0,R)
                gene=sample(seq(1,R,1),1)
                deltas[gene]=rnorm(1,mean=0,sd=Sigma) ## Normal distribution
                mutant=pmax(rep(0,R),parental+deltas)
              }
              # update genotypes counts
              Nt=Nt+1
              trackID=c(trackID,Nt)
              State=c(State,n0)
              t_alphas=rbind(t_alphas,mutant)
            }
            
             State[x]=max(0,State[x]-exp_mut) ## update densities
          }
            
        }
        tN=N+(Nt-Ns) ## temporary number of strains
        enzymes=unlist(lapply(1:R,function(x)(sum(State*t_alphas[,x])))) ## compute the state of the ecosystem
        growths=unlist(lapply(1:tN, function(s) sum(unlist(lapply(1:R,function(x) t_alphas[s,x]*sS[x]/enzymes[x]))))) # compute the expected growth
        exTout=pmax(rep(0,tN),State*(1+growths-d)) #expected update of the density
        
        ## multinomial sampling from the expected densities
        Tout=rmultinom(1,sum(exTout),exTout/sum(exTout)) ## Selection + Drift
        #Tout=rmultinom(1,sum(State),State/sum(State)) ## Pure drift
        #Tout=exTout ## Pure selection
        
        densities=as.numeric(Tout)
        surv_id=which(densities>=n0) ##eliminate the extinct
        New=setdiff(surv_id,1:N) ## get the new genotypes
        old=setdiff(surv_id,New)
          
        ## update the philogenetyc info
        trackID=trackID[old]
        if(length(New)>0){
          for(new in 1:length(New)){
            trackID=c(trackID,Ns+new)
            parent=parents[New[new]-N]
            history[[Ns+new]]=c(Ns+new,history[[parent]])
          }
        }
        ## update overall alphas
        alphas=rbind(alphas,t_alphas[New,])
        ## update genotypes count and densities
        N=length(surv_id)
        Ndym=c(Ndym,N)
        Ns=Ns+length(New)
        State=densities[surv_id]
        Dout[[generation]]=rbind(trackID,State)
          
        ## Compute Fitness proxy
        av_alphas=c(av_alphas,sum(State*rowSums(t_alphas[surv_id,]))/sum(State))
      }
      ## update data matrices
      Strains=rbind(get0("Strains"),Ndym)
      Av_alphas=rbind(get0("Av_alphas"),av_alphas)
        
      Alphas[[replica]]=alphas
      Dyn[[replica]]=Dout
      History[[replica]]=history
      Surv[[replica]]=trackID
      rm(list=c("alphas","Dout","history","trackID"))
    }
      
  ## compute mean and sd across simulations
  for (name in My_data){
    ## mean over simulations
    assign(name,rbind(get(name),colMeans(get(name))))
    ## 2*standard error
    assign(name,rbind(get(name),apply(get(name)[1:simulations,],2,sd)/sqrt(simulations)*2))
  }
  
  CPtime_table[1,u]=(proc.time()[3]-pt0)/simulations      
  print(paste("It took:",round(CPtime_table[1,u],digits=3), " seconds" , "(= ",round(CPtime_table[1,u]/3600,digits=3)," hours)  per simulation."))
  print("#####################################################################")
      
  allStrains[[counter]]=Strains
  allAv_alphas[[counter]]=Av_alphas
      
  rm(list=My_data)
  allAlphas[[counter]]=Alphas
  allDyn[[counter]]=Dyn
  allHistory[[counter]]=History
  allSurv[[counter]]=Surv
  rm(list=c("Alphas","Dyn","History","Surv"))
  ## save every iteration
  save.image(file=OutFile)
}
 
allCPtime=CPtime_table
save.image(file=OutFile)

