
#library(Rfast)

library(readr)
library(reshape2)
library(Rcpp)
library(vegan)
library(parallel)
library(ggplot2)
#library(tidyverse)
calc_rich<-function(population){
  population[population>0]<-T
  colSums(population)
}
fast_mean <- function(x) {
  x<-x[!is.na(x)]
  sum(x) / length(x)
} 

calc_div<-function(population){
  relabund<-population/rowSums(population)
  - rowSums(relabund* log(relabund))
}

fitness_func <- function(selection_parameter, optima, trait) {
  q = exp(((trait - optima)^2) / -selection_parameter)
  q
}


fitness_func_bacgen <- function(selection_parameter,optima1, optima2, trait) {
  q = exp((((trait - mean(c(optima1, optima2)))^2))/ -selection_parameter)
  q
}

mutate_trait <- function(trait, mutation_rate, mutation_sd) {
  mutated_trait <- trait #+ rnorm(length(trait), mean = 0, sd = mutation_sd)
  mutate_mask <- runif(length(trait)) < mutation_rate
  mutated_trait[mutate_mask] <- rnorm(sum(mutate_mask), mean = mutated_trait[mutate_mask], sd = mutation_sd)
  mutated_trait
}


Rcpp::sourceCpp('./calccpp_fit.cpp')
sourceCpp("./withingen_process.cpp")
sourceCpp("./env_production_withingen.cpp")

Rcpp::sourceCpp("./create_offspring_loop.cpp")

Create_OffSpring_Pop_cpp<-function(host_pop,n_micro,env_pool,envpoolsize,X,fixed_envpool,
                                   selection_parameter_microbes,selection_parameter_hosts,microbe_trait_list, host_microbe_optima,
                                   N_Species,env_condition){
  
  microbe_names<-paste("Microbe",1:N_Species,sep="_")
  #offspring_population<-list()
  #  fitness_microbes<-matrix(data=NA,nrow=N_Species,ncol=ncol(host_pop))
  ENV_sampling_probability<-(env_pool)/envpoolsize #Calculate initial sampling probability based on environmental relative abundance
  #  fitness_microbes<-vector()
  #  weighted_samplingprob<-vector()
  names(ENV_sampling_probability)<-names(env_pool)
  offspring_population<-offspring_loopfunc(host_pop = host_pop,ENV_sampling_probability = ENV_sampling_probability,
                                           host_microbe_optima = host_microbe_optima, X = X,
                                           selection_parameter_microbes = selection_parameter_microbes,
                                           env_condition = env_condition,microbe_trait_list = microbe_trait_list,
                                           n_micro = n_micro,microbe_names = microbe_names)#,
  #                                          fitness_microbes = fitness_microbes,weighted_samplingprob = weighted_samplingprob)
  offspring_population<-list(offspring_population$host_pop,env_pool,microbe_trait_list,offspring_population$fitness_microbes,
                             offspring_population$weighted_samplingprob)
  
  names(offspring_population)<-c("Child","Env","microbe_trait_list","microbefitness","microbe_samplingprob")
  offspring_population
}


lapply_wrapper_CPP<-function(XY,
                             HostPopulation, N_Microbes, envpoolsize, env_pool,fixed_envpool,generations,per_host_bac_gens,self_seed_prop,
                             selection_parameter_hosts,selection_parameter_microbes,host_trait_list,N_Species,traitpool_microbes,
                             generation_data_file,selection_parameter_env,env_cond_val,
                             microbiome_importances,host_microbe_optima,mutation_rate,mutation_sd,print_currentgen){
  
  # profvis({
  #Note: selection_parameter_hosts can be used to turn on, and adjust intensity of HS - when zero, host selection is off.
  #Similarly when selection_parameter_microbes is 0, we are using a neutral model, regardless of the given value of SelectionType. 
  temp_list <- list()
  print(paste("Current X Value is", XY[1], "Current EnvCon value is", XY[2]))
  gen_data<-list()    
  microbe_names<-paste("Microbe",1:N_Species,sep="_")
  #print(dim(env_cond_val))
  #print((env_cond_val))
  
  env_used <- NULL
  if(!nrow(env_cond_val) == generations){
    print("env_cond_val should be a matrix of the same length as the number of generations you are simulating")
    stop()
  }
  init_microbial_fitness<-matrix(NA,nrow = length(traitpool_microbes),ncol=length(host_microbe_optima))
  
  for(i in 1:ncol(init_microbial_fitness)){
    init_microbial_fitness[,i] <- fitness_func_bacgen(selection_parameter = selection_parameter_microbes,
                                                      optima1 = host_microbe_optima[i],
                                                      optima2 = env_cond_val[1,1],trait = traitpool_microbes)
  }
  init_microbial_fitness<-init_microbial_fitness[HostPopulation>0]
  
  #  relabund<-t(apply(HostPopulation,1,function(x){x/N_Microbes}))
  init_microbial_fitness<-weighted.mean(x = init_microbial_fitness,w = HostPopulation)
  
  
  init_host_fitness<-vector()
  for(i in 1:ncol(HostPopulation)){
    mean_microbial_trait_val<- sum(traitpool_microbes * HostPopulation[,i])/N_Microbes
    composite_host_trait<-((mean_microbial_trait_val * microbiome_importances[i]) + (host_microbe_optima[i] * (1-microbiome_importances[i])))
    hostfitness<-fitness_func(selection_parameter = selection_parameter_hosts, trait = composite_host_trait,optima = env_cond_val[1,1])
    init_host_fitness[i]<-hostfitness
  }
  
  if(!ncol(env_cond_val)==per_host_bac_gens){
    print("Number of columns in env_conditions should be the number of bacterial generations")
    stop()
  }
  if(is.na(self_seed_prop)){
    print("Self seeding proportion is NA, please correct")
  }
  for (generation in 1:generations) {
    if(generation ==1){
      env_used<-fixed_envpool
    }
#    if(per_host_bac_gens == 1){
#      self_seed_prop=1
#    }
    env_cond_val_used<-as.vector(env_cond_val[generation,])
    # print(env_cond_val_used)
    #  print(paste("Current Generation is", generation))
    for(bacgen in 1:per_host_bac_gens){
      
      if(XY[2] + XY[3] > 1) {
        print("The contribution of the fixed environment (Y) and the the autocthonous environment (var_env_con) is greater than 1, please correct this ")
        stop()
      }
      #print(min(env_used))
      gen_env_cond<-env_cond_val_used[bacgen]
      #   print(gen_env_cond)
      env_used<-process_microbe_probs(env_cond_val = gen_env_cond,
                                      fixed_envpool = fixed_envpool,
                                      HostPopulation = HostPopulation,
                                      N_Microbes = N_Microbes,
                                      envpoolsize = envpoolsize,
                                      selection_parameter_env = selection_parameter_env,
                                      XY = XY,
                                      traitpool_microbes = traitpool_microbes,env_used = env_used
      )#[,8]
      env_fits<-weighted.mean(x = env_used[,3],w = env_used[,8])
      env_used<-env_used[,8]        
      names(env_used)<-names(fixed_envpool)#rownames(microbe_probs)
      # print(table(is.na(env_used)))
        if(bacgen >1){
      HostPopulation<-process_host_cpp(HostPopulation = HostPopulation,
                                       selection_parameter_microbes = selection_parameter_microbes,
                                       host_microbe_optima = host_microbe_optima ,
                                       env_condition = env_cond_val_used[bacgen],
                                       traitpool_microbes = traitpool_microbes,
                                       N_Microbes = N_Microbes, 
                                       self_seed_prop = self_seed_prop,env_used = env_used)}
        else {
        HostPopulation<-HostPopulation
            }
      #result <- lapply(seq_along(colnames(HostPopulation)), process_host)
      #print(HostPopulation)
    }
    #Now we've created our new environments, we now need to choose which members of a population reproduce
    # We can do this either neutrally (i.e., random chance) OR we can do this based on the fitness of a host based on what is provided by its microbiome
    host_fitnessvector<-numeric()
    #  print(range(HostPopulation))
    host_fitnessvector <- calculate_host_fitness_cpp(as.matrix(HostPopulation), 
                                                     traitpool_microbes, microbiome_importances, 
                                                     host_microbe_optima, 
                                                     selection_parameter_hosts, env_cond_val_used[per_host_bac_gens])
    nomicrobiomehost_fitnessvector<-numeric()
    nomicrobiomehost_fitnessvector <- calculate_host_fitness_cpp(as.matrix(HostPopulation), 
                                                                 traitpool_microbes,rep(0,length(microbiome_importances)), 
                                                                 host_microbe_optima, 
                                                                 selection_parameter_hosts, env_cond_val_used[per_host_bac_gens])
    
    #    names(host_fitnessvector)<-colnames(HostPopulation)
    #    
    hostfitness_abs<-host_fitnessvector
    # print(host_fitnessvector)
    host_fitnessvector<-host_fitnessvector#/mean(host_fitnessvector)
    
    HostPopulationInt<-sample(colnames(HostPopulation),ncol(HostPopulation),replace = T,prob = host_fitnessvector)
    #host_preference_optima_prevgen<-host_preference_optima
    host_microbe_optima_prevgen<-host_microbe_optima
    
    HostPopulation<-HostPopulation[ , HostPopulationInt]
    host_microbe_optima<-host_microbe_optima[HostPopulationInt]
    microbiome_importances<-microbiome_importances[HostPopulationInt]
    mean_microbe_traitval<-mean(apply(HostPopulation,MARGIN = 2,FUN=function(x){weighted.mean(traitpool_microbes,x)}))
    # host_preference_optima<-host_preference_optima[HostPopulationInt]
    
    if(mutation_rate>0){
      #  host_preference_optima <- mutate_trait(host_preference_optima, mutation_rate, mutation_sd)
      host_microbe_optima <- mutate_trait(host_microbe_optima, mutation_rate, mutation_sd)
      #  microbiome_importances <-
    }
    
    
    colnames(HostPopulation)<-paste("Host",1:ncol(HostPopulation),sep="_")
    names(host_microbe_optima)<-paste("Host",1:ncol(HostPopulation),sep="_")
    names(microbiome_importances)<-paste("Host",1:ncol(HostPopulation),sep="_")
    #  names(host_preference_optima)<-paste("Host",1:ncol(HostPopulation),sep="_")
    if(print_currentgen==T){print(paste("Current Generation is", generation))}
    new_gen <- Create_OffSpring_Pop_cpp(host_pop = HostPopulation,
                                        n_micro = N_Microbes,
                                        env_pool = env_used,
                                        envpoolsize = envpoolsize,
                                        X = XY[1],
                                        fixed_envpool = fixed_envpool, 
                                        microbe_trait_list=traitpool_microbes,
                                        selection_parameter_microbes=selection_parameter_microbes,
                                        selection_parameter_hosts = selection_parameter_hosts,
                                        host_microbe_optima=host_microbe_optima,
                                        N_Species=N_Species,
                                        env_condition=env_cond_val_used[per_host_bac_gens])
    # print("Issue here")
    HostPopulation <- new_gen$Child
    BrayDiv<-vegdist(t(HostPopulation),method="bray")
    div_pergen<-calc_div(t(new_gen$Child))/log(calc_rich(new_gen$Child))
    new_gen$HostFitness_Abs<-(hostfitness_abs)
    new_gen$nomicrobiomehost_fitnessvector<-nomicrobiomehost_fitnessvector
    new_gen$BrayDiv<-BrayDiv
    # print(new_gen$microbefitness[!is.na(new_gen$microbefitness)])
    #  new_gen$MicrobeFitness<-(new_gen$microbefitness)#,na.rm=T)#This values is the ACTUAL mean 
    #absolute fitness of the community at this time point. We earlier exclude absent communitiy members
    #new_gen$HostPrefOptima<-host_preference_optima_prevgen
    new_gen$HostMicrobeOptima<-host_microbe_optima_prevgen
    #is_coalesc <- identical(new_gen$HostMicrobeOptima, new_gen$HostMicrobeOptima)
    # print(new_gen$microbefitness)
    gen_data[[generation]]<-list(div_pergen,new_gen$HostFitness_Abs,new_gen$microbefitness,env_used,
                                 env_fits,new_gen$HostMicrobeOptima,new_gen$nomicrobiomehost_fitnessvector,new_gen$microbe_samplingprob,
                                 microbiome_importances,new_gen$BrayDiv,mean_microbe_traitval)
    names(gen_data[[generation]])<-c("Diversity","HostFitness","MicrobeFitness","env_used","env_fits"
                                     ,"HostMicrobeOptima","nomicrobiomehost_fitnessvector",
                                     "microbe_samplingprob","microbiome_importances","BrayDiv","mean_microbe_traitval")#,"HostPrefOptima","HostMicrobeOptima")
  }
  new_gen$GenData<-gen_data
  new_gen$init_host_fitness<-init_host_fitness
  new_gen$init_microbe_fitness<-init_microbial_fitness
  new_gen$traitpool<-traitpool_microbes
  #  if(!is.na(generation_data_file)){
  #    readr::write_rds(gen_data,generation_data_file,compress = "gz")
  #  }
  temp_list[[paste0("X", XY[1], "_Y", XY[2],"_EnvCont",XY[3],"_HostSel",XY[4],"_MicrobeSel",XY[5])]] <- new_gen
  temp_list
  #})
}
num_unique_values <- function(vec) {
  length(unique(vec))
}

# Function to calculate average distance between values
average_distance <- function(vec) {
  n <- length(vec)
  total_distance <- 0
  total_pairs <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      total_distance <- total_distance + abs(vec[i] - vec[j])
      total_pairs <- total_pairs + 1
    }
  }
  total_distance / total_pairs
}
extract_pref_optima<-function(x){
  do.call("rbind",lapply(x[[1]]$GenData,function(y){(y$HostPrefOptima)}))
}
extract_fitnesses<-function(data,rel_log_gen){
  gendata=data$GenData
  data2use=lapply(gendata,function(q){
    env_microbe_fitness_changes<-mean((q$env_fits))#/x[[1]]$GenData[[500]]$env_fits))
    host_microbe_fitness_changes<-mean((q$MicrobeFitness))#/x[[1]]$GenData[[500]]$MicrobeFitness))
    host_fitness_changes<-mean((q$HostFitness))#/x[[1]]$GenData[[500]]$HostFitness))
    gen_of_int_fitness_host<-mean(gendata[[rel_log_gen]]$HostFitness)
    gen_of_int_fitness_microbe<-mean(gendata[[rel_log_gen]]$MicrobeFitness)
    host_microbe_coeff=lm(q$MicrobeFitness~q$HostFitness)
    host_microbe_coeff<-summary(host_microbe_coeff)$coefficients
    host_microbe_coeff<-as.data.frame(host_microbe_coeff)$Estimate[2]
    host_diversity<-mean(q$Diversity)
    data.frame(EnvMicrobeFitness=env_microbe_fitness_changes,
               HostMicrobeFitness=host_microbe_fitness_changes,
               HostFitness=host_fitness_changes,
               diversity=host_diversity,
               gen_of_int_fitness_microbe=gen_of_int_fitness_microbe,
               gen_of_int_fitness_host=gen_of_int_fitness_host,
               host_microbe_coeff=host_microbe_coeff)
  }
  )
  do.call("rbind",data2use)
}


##############

rep_number <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#for(rep_number in 2:20){
Host_PopSize= 100
MicrobePopSize=10^9#0000 #n_micro
EnvPoolSize=10^9#0000
N_Species =200
N_Traits=25
N_HostTraits=25
N_Microbe_Traits=5
selection_type="HMS"
InitPopulation<-matrix(nrow = N_Species, ncol=Host_PopSize)
InitPopulation<-apply(InitPopulation,MARGIN=2,FUN=function(x){
  rmultinom(n=1,size = MicrobePopSize, prob =rep(1/N_Species,N_Species))
})
rownames(InitPopulation)<-paste("Microbe",1:N_Species,sep="_")
colnames(InitPopulation)<-paste("Host",1:Host_PopSize,sep="_")

EnvPool<-rep(EnvPoolSize/N_Species,N_Species)
names(EnvPool)<-paste("Microbe",1:N_Species,sep="_")
fixed_envpool<-EnvPool
rnames<-paste("Host",1:Host_PopSize,sep="_")
cnames<-paste("Trait",1:N_HostTraits,sep="_")
rnames<-paste("Microbe",1:N_Species,sep="_")
cnames<-paste("Trait",1:N_Traits,sep="_")
traitpool_microbes<-runif(N_Species,-1,1)#truncnorm::rtruncnorm(n = 200,-1,1,mean = 0,sd = 0.6)#runif(200,-1,1)##rowMeans(traitpool_microbes,na.rm=T)
names(traitpool_microbes)<-paste("Microbe",1:N_Species,sep="_")
mnames<-colnames(InitPopulation)
host_microbe_optima<-rep(0,Host_PopSize)##truncnorm::rtruncnorm(n = Host_PopSize,-1,1,mean = 0,sd = 0)
#host_microbe_optima<-rep(runif(1,-1,1),Host_PopSize)##truncnorm::rtruncnorm(n = Host_PopSize,-1,1,mean = 0,sd = 0)


names(host_microbe_optima)<-mnames
microbiome_importances<-rep(1,Host_PopSize)#runif(Host_PopSize,0,1)
names(microbiome_importances)<-mnames
nhost_gens<-1501
N=1501
burn_length <- 200#round(N * burn_percentage)
burn_noise <- rnorm(burn_length, mean =0, sd = 0.1)
transition_length <- 1100
transition_values <- seq(0, 0.6, length.out = transition_length)
transition_noise <- rnorm(transition_length, mean = 0, sd = 0.1)
remaining_length <- N - burn_length - transition_length
remaining_noise <- rnorm(remaining_length, mean = 0.6, sd = 0.1)
vector1 <- c(burn_noise, transition_values + transition_noise, remaining_noise)
vector2 <- rnorm(N, mean = 0, sd = 0.1)
N=nhost_gens
start <- seq(0.1, 0.1, length.out = N*0.2)
mid <- seq(0.1, 0.3, length.out = N*0.8)
end <- seq(0.3, 0.3, length.out = N*0.2)
library(stats)
simulate_arima_varying_sd <- function(ar, ma, n, sd_vector) {
  if(length(sd_vector) != n) {
    stop("Length of sd_vector must be equal to n")
  }
  #
  #   # Initialize the time series
  time_series <- numeric(n)
  #
  #   # Generate white noise with varying standard deviations
  white_noise <- rnorm(n, mean = 0, sd = sd_vector)
  #
  #   # Simulate the ARIMA process manually
  for (t in 1:n) {
    if (t == 1) {
      time_series[t] <- white_noise[t]
    } else {
      ar_part <- ifelse(length(ar) > 0, sum(ar * head(time_series, t)[(t-length(ar)):(t-1)]), 0)
      ma_part <- ifelse(length(ma) > 0, sum(ma * head(white_noise, t)[(t-length(ma)):(t-1)]), 0)
      time_series[t] <- ar_part + ma_part + white_noise[t]
    }
  }
  #
  return(time_series)
}
highac_incvar <- simulate_arima_varying_sd(0.9, 0, nhost_gens, c(start,mid)[1:1501])
lowac_incvar<-sapply(c(start,mid)*2.8, function(x){rnorm(1,0,x)})
change_env_vec<-c(rep(0,length=0.2*(nhost_gens+1)),seq(0,1,length.out=0.8*(nhost_gens+1)))
highac <- simulate_arima_varying_sd(0.9, 0, (nhost_gens+1),rep(0.2,(nhost_gens+1)))#
lowac<-sapply(rep(0.2,nhost_gens+1)*2.6, function(x){rnorm(1,0,x)})
highac_incmean <- simulate_arima_varying_sd(0.9, 0, (nhost_gens+1), rep(0.2,(nhost_gens+1)))+change_env_vec
lowac_incmean <-sapply(rep(0.2,(nhost_gens+1))*2.3, function(x){rnorm(1,0,x)})+change_env_vec
n <- nhost_gens+1  # Total length of the vector including burn-in
burn_in <- 200  # Length of burn-in period
transition_length <- n - burn_in  # Length of transition period
autocorr_values_burn_in <- rep(0.5, burn_in)  # Autocorrelation values during burn-in
autocorr_values_transition <- seq(0.5, 0.99, length.out = transition_length)  # Autocorrelation values during transition
autocorr_values <- c(autocorr_values_burn_in, autocorr_values_transition)  # Combined autocorrelation values
sd <- 0.2  # Standard deviation
interpolate_vector <- function(vec, N) {
  # Calculate total number of points to interpolate
  new_length <- (length(vec) - 1) * (N + 1) + 1
  
  # Perform interpolation
  interpolated_data <- approx(x = 1:length(vec), y = vec, n = new_length)$y
  
  return(interpolated_data)
}
expand_vector<-function(nhostgen,nmicrogen,env_vec){
  newvec<-interpolate_vector(env_vec,nmicrogen-1)
  newvec<-newvec[1:(nmicrogen*nhostgen)]
  newvec<-matrix(newvec,nrow=nhostgen,ncol=nmicrogen,byrow=T)
  newvec
  #table(newvec[,1]==lowac[1:nhostgen])
}
lowac_expand<-list()
lowac_incmean_expand<-list()
lowac_incvar_expand<-list()
highac_expand<-list()
highac_incmean_expand<-list()
highac_incvar_expand<-list()
for(i in c(1,2,3,4,5,6,7,8,9,10,20,50,100,200)){
  lowac_expand[[i]]<-expand_vector(1500,i,lowac)
  lowac_incmean_expand[[i]]<-expand_vector(1500,i,lowac_incmean)
  lowac_incvar_expand[[i]]<-expand_vector(1500,i,lowac_incvar)
  highac_expand[[i]]<-expand_vector(1500,i,highac)
  highac_incmean_expand[[i]]<-expand_vector(1500,i,highac_incmean)
  highac_incvar_expand[[i]]<-expand_vector(1500,i,highac_incvar)
}
varying_envs<-list(lowac_expand,lowac_incmean_expand,lowac_incvar_expand,highac_expand,highac_incmean_expand,highac_incvar_expand)
names(varying_envs)<-c("lowac","lowac_incmean","lowac_incvar","highac","highac_incmean","highac_incvar")
varying_envs<-lapply(varying_envs,function(x){
  x<-x[!unlist(lapply(x,is.null))]
  x
})
varying_envs<-lapply(varying_envs,function(x){
  names(x)<-paste0(c(1,2,3,4,5,6,7,8,9,10,20,50,100,200),"gen")
  x
})


VerticalInheritance<-c(0,16.67,33.34,50,66.7,93.35,99,99.9,100)/100

HostCont= c(5)/100
EnvSelf = c(80)/100
HostSelStrength=c(1)
MicrobeSelStrength=c(1)
EnvSelStrength=c(1)
sec_varyring_envs_multigen<-unlist(varying_envs,recursive = F)

write_rds(x = sec_varyring_envs_multigen,paste0("rerun_1nov_rerun_sec_varyring_envs_multigen","_replicate_",rep_number,".RDS"),compress = "gz")


all_combinations<-expand.grid(VerticalInheritance,EnvSelf,HostCont,HostSelStrength,MicrobeSelStrength,EnvSelStrength,names(sec_varyring_envs_multigen))
run_simulation <- function(combination) {
  env_matrix <- eval(parse(text=paste0("sec_varyring_envs_multigen$",as.character(combination$Var7))))
  combination<-as.vector(unlist(combination))
  combination<-as.numeric(combination[1:6])
  env_matrix<-as.matrix(env_matrix)
  lapply_wrapper_CPP(
    XY = combination,
    HostPopulation = InitPopulation,
    N_Microbes = MicrobePopSize,
    envpoolsize = EnvPoolSize,
    env_pool = EnvPool,
    generations = nrow(env_matrix),
    fixed_envpool = fixed_envpool,
    selection_parameter_hosts = combination[4],
    selection_parameter_microbes = combination[5],
    traitpool_microbes = traitpool_microbes,
    host_microbe_optima = host_microbe_optima,
    selection_parameter_env = combination[6],
    env_cond_val = env_matrix,
    N_Species = N_Species,
    microbiome_importances = microbiome_importances,
    per_host_bac_gens = ncol(env_matrix),
    self_seed_prop = 0.98,
    generation_data_file = NA,
    mutation_rate = 0.0,
    mutation_sd = 0.0,
    print_currentgen = F
  )
}
results_multigen <- mclapply(seq_len(nrow(all_combinations)), function(i) {
  run_simulation(all_combinations[i, ])
}, mc.cores = 60)
results_multigen_unlist<-unlist(results_multigen,recursive = F)
names(results_multigen_unlist)<-(paste(names(results_multigen_unlist),all_combinations$Var7,sep="_"))
write_rds(x=results_multigen_unlist,file = paste0("rerun_1nov_largererun_results_multigen_replicate_",rep_number,".RDS"),compress = "gz" )
rm(results_multigen_unlist)
rm(results_multigen)
