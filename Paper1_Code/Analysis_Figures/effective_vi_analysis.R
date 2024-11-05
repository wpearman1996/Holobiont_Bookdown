setwd("~/Dropbox/Anna_PostDoc/simulation_code_and_results/holobiont_paper1/evi_analysis//")

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


Rcpp::sourceCpp('calccpp_fit.cpp')
sourceCpp("withingen_process.cpp")
sourceCpp("env_production_withingen.cpp")

Rcpp::sourceCpp("create_offspring_loop_vi.cpp")

Create_OffSpring_Pop_cpp<-function(host_pop,n_micro,env_pool,envpoolsize,X,fixed_envpool,
                                   selection_parameter_microbes,selection_parameter_hosts,microbe_trait_list, host_microbe_optima,
                                   N_Species,env_condition,generation){
  
  microbe_names<-paste("Microbe",1:N_Species,sep="_")
  #offspring_population<-list()
  #  fitness_microbes<-matrix(data=NA,nrow=N_Species,ncol=ncol(host_pop))
  ENV_sampling_probability<-(env_pool)/envpoolsize #Calculate initial sampling probability based on environmental relative abundance
  #  fitness_microbes<-vector()
  #  weighted_samplingprob<-vector()
  names(ENV_sampling_probability)<-names(env_pool)
  if(generation == 498){
    offspring_population<-offspring_loopfunc_vi(host_pop = host_pop,ENV_sampling_probability = ENV_sampling_probability,
                                                host_microbe_optima = host_microbe_optima, X = X,
                                                selection_parameter_microbes = selection_parameter_microbes,
                                                env_condition = env_condition,microbe_trait_list = microbe_trait_list,
                                                n_micro = n_micro,microbe_names = microbe_names)
    # print(head(offspring_population$host_pop,100))
  }
  else{
    offspring_population<-offspring_loopfunc(host_pop = host_pop,ENV_sampling_probability = ENV_sampling_probability,
                                             host_microbe_optima = host_microbe_optima, X = X,
                                             selection_parameter_microbes = selection_parameter_microbes,
                                             env_condition = env_condition,microbe_trait_list = microbe_trait_list,
                                             n_micro = n_micro,microbe_names = microbe_names)#,
    #                                          fitness_microbes = fitness_microbes,weighted_samplingprob = weighted_samplingprob)
  }
  offspring_population<-list(offspring_population$host_pop,env_pool,microbe_trait_list,offspring_population$fitness_microbes,
                             offspring_population$weighted_samplingprob)
  
  names(offspring_population)<-c("Child","Env","microbe_trait_list","microbefitness","microbe_samplingprob")
  offspring_population
}


lapply_wrapper_CPP<-function(XY,
                             HostPopulation, N_Microbes, envpoolsize, env_pool,fixed_envpool,generations,per_host_bac_gens,self_seed_prop,
                             selection_parameter_hosts,selection_parameter_microbes,host_trait_list,N_Species,traitpool_microbes,
                             generation_data_file,selection_parameter_env,env_cond_val,
                             microbiome_importances,host_microbe_optima,mutation_rate,mutation_sd,print_currentgen,env_type){
  
  # profvis({
  #Note: selection_parameter_hosts can be used to turn on, and adjust intensity of HS - when zero, host selection is off.
  #Similarly when selection_parameter_microbes is 0, we are using a neutral model, regardless of the given value of SelectionType.
  temp_list <- list()
  print(paste("Current X Value is", XY[1], "Current EnvCon value is", XY[2]))
  print(env_type)
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
  init_microbial_fitness<-init_microbial_fitness#[HostPopulation>0]
  # HostPopulation<-as.matrix(HostPopulation)
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
    env_cond_val_used<-as.vector(env_cond_val[generation,])
    if(generation == 499){
      print(HostPopulation[,1])
    }
    #    HostPopulation[,1]<-as.vector(c(HostPopulation[201:400,1],HostPopulation[1:200,1]))
    vi_ests<-matrix(NA,nrow = 400,ncol = 200)#list()
    # }
    for(bacgen in 1:per_host_bac_gens){
      
      if(XY[2] + XY[3] > 1) {
        print("The contribution of the fixed environment (Y) and the the autocthonous environment (var_env_con) is greater than 1, please correct this ")
        stop()
      }
      #print(min(env_used))
      #print(env_used)
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
      # if(generation == 499){
      #   vi_ests[[bacgen]]<-HostPopulation
      # }
      if(bacgen > 1){
      HostPopulation<-process_host_cpp(HostPopulation = HostPopulation,
                                       selection_parameter_microbes = selection_parameter_microbes,
                                       host_microbe_optima = host_microbe_optima ,
                                       env_condition = env_cond_val_used[bacgen],
                                       traitpool_microbes = traitpool_microbes,
                                       N_Microbes = N_Microbes,
                                       self_seed_prop = self_seed_prop,env_used = env_used)
      }
      
      if(generation == 499){
        #print(max(HostPopulation[,1]))
        vi_ests[,bacgen]<-HostPopulation[,1]
      }
      #print(HostPopulation)
      #result <- lapply(seq_along(colnames(HostPopulation)), process_host)
      #   print(bacgen)
    }
    if(generation == 499){
      write_rds(vi_ests,paste("X",env_type,XY[1],"vi_ests.RDS",sep="_"))
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
    
    # host_preference_optima<-host_preference_optima[HostPopulationInt]
    
    if(mutation_rate>0){
      #  host_preference_optima <- mutate_trait(host_preference_optima, mutation_rate, mutation_sd)
      host_microbe_optima <- mutate_trait(host_microbe_optima, mutation_rate, mutation_sd)
      #  microbiome_importances <-
    }
    
    HostPopulation<-as.matrix(HostPopulation)
    # if(generation == 499){
    #   print(vi_ests)
    # }
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
                                        env_condition=env_cond_val_used[per_host_bac_gens],
                                        generation = generation)
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
                                 microbiome_importances,new_gen$BrayDiv)
    names(gen_data[[generation]])<-c("Diversity","HostFitness","MicrobeFitness","env_used","env_fits"
                                     ,"HostMicrobeOptima","nomicrobiomehost_fitnessvector",
                                     "microbe_samplingprob","microbiome_importances","BrayDiv")#,"HostPrefOptima","HostMicrobeOptima")
  }
  new_gen$GenData<-gen_data
  # new_gen$vi_ests<-vi_ests
  new_gen$init_host_fitness<-init_host_fitness
  new_gen$init_microbe_fitness<-init_microbial_fitness
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
Host_PopSize= 100
MicrobePopSize=10^9#0000 #n_micro
EnvPoolSize=10^9#0000
N_Species =400
N_Traits=25
N_HostTraits=25
N_Microbe_Traits=5
selection_type="HMS"
InitPopulation<-matrix(nrow = 200, ncol=Host_PopSize)
InitPopulation<-apply(InitPopulation,MARGIN=2,FUN=function(x){
  rmultinom(n=1,size = MicrobePopSize, prob =rep(1/200,200))
})
rownames(InitPopulation)<-paste("Microbe",1:200,sep="_")
colnames(InitPopulation)<-paste("Host",1:Host_PopSize,sep="_")


VI_InitPopulation<-matrix(data=0,nrow = 200, ncol=Host_PopSize)
rownames(VI_InitPopulation)<-paste("VI_Microbe",1:200,sep="_")
colnames(VI_InitPopulation)<-paste("Host",1:Host_PopSize,sep="_")
InitPopulation<-rbind(InitPopulation,VI_InitPopulation)


EnvPool<-rep(EnvPoolSize/200,200)
names(EnvPool)<-paste("Microbe",1:200,sep="_")

VI_EnvPool<-rep(0,200)
names(VI_EnvPool)<-paste("VI_Microbe",1:200,sep="_")
EnvPool<-c(EnvPool,VI_EnvPool)

fixed_envpool<-EnvPool




rnames<-paste("Host",1:Host_PopSize,sep="_")
cnames<-paste("Trait",1:N_HostTraits,sep="_")
rnames<-paste("Microbe",1:N_Species,sep="_")
cnames<-paste("Trait",1:N_Traits,sep="_")
traitpool_microbes<-runif(200,-1,1)#truncnorm::rtruncnorm(n = 200,-1,1,mean = 0,sd = 0.6)#runif(200,-1,1)##rowMeans(traitpool_microbes,na.rm=T)
traitpool_microbes<-c(traitpool_microbes,traitpool_microbes)
names(traitpool_microbes)<-names(EnvPool)
mnames<-colnames(InitPopulation)
#host_microbe_optima<-rep(rnorm(n = 1,0,0.2),Host_PopSize)##truncnorm::rtruncnorm(n = Host_PopSize,-1,1,mean = 0,sd = 0)
host_microbe_optima<-rep(0,Host_PopSize)##truncnorm::rtruncnorm(n = Host_PopSize,-1,1,mean = 0,sd = 0)



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
sec_varyring_envs_multigen<-unlist(varying_envs,recursive = F)
sec_varyring_envs_multigen<-sec_varyring_envs_multigen[grepl("200gen",names(sec_varyring_envs_multigen))]
env_new_vectors<-lapply(sec_varyring_envs_multigen,function(x){x[1:500,]})

VerticalInheritance<-c(0,16.67,33.34,50,66.7,93.35,99,99.9,100)/100
HostCont= c(5)/100
EnvSelf = c(80)/100
HostSelStrength=c(1)
MicrobeSelStrength=c(1)
EnvSelStrength=c(1)
#sec_varyring_envs_multigen<-env_new_vectors[c("lowac","highac_incvar","highac_incmean","lowac","lowac_incvar","lowac_incmean")]
#sec_varyring_envs_multigen<-env_new_vectors[c("lowac")]
#sec_varyring_envs_multigen<-env_new_vectors[c("lowac_incvar","lowac_incmean","highac")] 

#sec_varyring_envs_multigen<-unlist(sec_varyring_envs_multigen,recursive = F)
#write_rds(x = sec_varyring_envs_multigen,paste0("altdist_wn_largererun_sec_varyring_envs_multigen","_replicate_",rep_number,".RDS"),compress = "gz")


all_combinations<-expand.grid(VerticalInheritance,EnvSelf,HostCont,HostSelStrength,MicrobeSelStrength,EnvSelStrength,names(sec_varyring_envs_multigen))
#all_combinations<-all_combinations[c(1:2,8:9),]
run_simulation <- function(combination) {
  env_matrix <- eval(parse(text=paste0("sec_varyring_envs_multigen$",as.character(combination$Var7))))
  envtype<-as.character(combination$Var7)
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
    print_currentgen = F,
    env_type=envtype
  )
}
results_multigen <- mclapply(seq_len(nrow(all_combinations)), function(i) {
  run_simulation(all_combinations[i, ])
}, mc.cores = 6)
#
results_multigen_unlist<-unlist(results_multigen,recursive = F)
names(results_multigen_unlist)<-(paste(names(results_multigen_unlist),all_combinations$Var7,sep="_"))
write_rds(x=results_multigen_unlist,file = paste0("effective_VI_analysis_4Nov",".RDS"),compress = "gz" )
# rm(results_multigen_unlist)
# rm(results_multigen)
viests<-list.files(pattern="X*ests.RDS")
viests<-lapply(viests,read_rds)
viests<-lapply(viests,function(q){
  colSums(q[201:400,])/colSums(q)
})

names(viests)<-list.files(pattern="X*ests.RDS")

viests<-do.call("rbind",viests)
viests<-reshape2::melt(viests)
viests<-viests[complete.cases(viests),]
viests$InitVert<- gsub(".RDS","",viests$Var1)
viests$InitVert<- stringr::word(viests$InitVert,2,sep="200gen")
viests$InitVert<-gsub("[^0-9.-]", "", viests$InitVert)

viests$Environment<-gsub("[[:digit:]]+","",viests$Var1)
viests$Environment<-gsub("._vi_ests.RDS","",viests$Environment)
viests$Environment<-gsub("X_","",viests$Environment)
viests$Environment<-gsub(".gen_","",viests$Environment)
efi_relationship<-ggplot(viests) +
  aes(x = Var2, y = value, colour = InitVert) +
  geom_line(linewidth=2) +
  scale_color_hue(direction = 1) +
  theme_minimal() + xlab("Microbial Generations within a host") + 
  ylab("Proportion of Vertical Acquired Microbes remaining") #+ facet_wrap(~Environment)
efi_relationship
plot(viests$value~as.numeric(viests$Var2),col=as.factor(viests$Var1))
viests[viests$Var2==200,]


all_combinations<-expand.grid(VerticalInheritance,EnvSelf,HostCont,HostSelStrength,MicrobeSelStrength,EnvSelStrength,names(sec_varyring_envs_multigen))
all_combinations_evi<-all_combinations
viests_int<-viests[viests$Var2==200,]
viests_int$Lookupval<-paste(viests_int$Environment,".200gen",viests_int$InitVert ,sep="")
all_combinations_evi$Lookup<-gsub("_","",all_combinations_evi$Var7)
all_combinations_evi$Lookup<-paste(all_combinations_evi$Lookup,all_combinations_evi$Var1,sep="")
all_combinations_evi$Var1<-viests_int$value[match(all_combinations_evi$Lookup,viests_int$Lookupval)]
all_combinations_evi$Var1<-viests$value[viests$Var2==200]
all_combinations_evi$Lookup<-NULL
all_combinations_evi$Var7<-gsub("200gen","1gen",all_combinations_evi$Var7)
#write.table(x = all_combinations_evi,"../evi_analysis/all_combinations_evi4nov.txt")

all_combinations_evi<-read.table("../../evi_analysis//all_combinations_evi4nov.txt")
all_combinations__orig<-all_combinations
all_combinations_evi$VertConOrig<-all_combinations__orig$Var1

######
#Effective sims analysis
library(readr)
library(ggpubr)
setwd("~/Dropbox/Anna_PostDoc/simulation_code_and_results/holobiont_paper1/draft_ms/evi_1gen_wn_largererun_data_analysed_replicate_4nov/")
rds_files<-list.files(pattern="*RDS")
rds_files<-lapply(rds_files,read_rds)
for(i in 1:length(rds_files)){
  rds_files[[i]]$Replicate<-paste0("replicate_",i)
}
effect_VI_data<-do.call("rbind",rds_files)
#levels_VerticalInheritance<-c(0,0.4636128,0.9240472, #27Aug - 200 gen update
#                       1.3777093,1.8308752,2.5637507,
#                       2.7190888,2.7423164,2.7442086)/100
#levels_VerticalInheritance<-as.character(levels_VerticalInheritance)

effect_VI_data$MeanGroupingColumn<-paste(effect_VI_data$Replicate,effect_VI_data$EnvType,
                                        paste0("1gen",effect_VI_data$Gen),sep=".")


effect_VI_data$MeanGroupingColumn<-paste(effect_VI_data$Gen,effect_VI_data$VertCon,
                                        effect_VI_data$EnvType,effect_VI_data$NGens)
effect_VI_data<-effect_VI_data[effect_VI_data$VertCon %in% all_combinations_evi$Var1,]


library(dplyr)
effect_VI_data_gen1<-effect_VI_data %>% filter(NGens == "1gen")
effect_VI_data_gen1$NewGrouping_Column<-gsub(" 1gen","",effect_VI_data_gen1$MeanGroupingColumn)
effect_VI_data_gen1$NewGrouping_Column<-paste(effect_VI_data_gen1$NewGrouping_Column,effect_VI_data_gen1$Replicate)
effect_VI_data$cleaned_strings <- gsub("\\S*gen\\S*", "", effect_VI_data$MeanGroupingColumn)
effect_VI_data$cleaned_strings <- trimws(effect_VI_data$cleaned_strings)
effect_VI_data$cleaned_strings<-paste(effect_VI_data$cleaned_strings,effect_VI_data$Replicate)

effect_VI_data$gen1_microbefitness<-effect_VI_data_gen1$HostMicrobeFitness[match(effect_VI_data$cleaned_strings,effect_VI_data_gen1$NewGrouping_Column)]
effect_VI_data$gen1_hostfitness<-effect_VI_data_gen1$HostFitness[match(effect_VI_data$cleaned_strings,effect_VI_data_gen1$NewGrouping_Column)]
effect_VI_data<-effect_VI_data %>% mutate(NGens=factor(NGens,levels=c("1gen","2gen","5gen","10gen","20gen","50gen","100gen","200gen")))

effect_VI_data$LogFold_ngen1_Host<-(effect_VI_data$HostFitness/effect_VI_data$gen1_hostfitness)
effect_VI_data$LogFold_ngen1_Microbe<-(effect_VI_data$HostMicrobeFitness/effect_VI_data$gen1_microbefitness)

effect_VI_data_means <- effect_VI_data %>%
   #filter(EnvType =="lowac") %>%
  group_by(EnvType,VertCon,NGens,Gen) %>%
  summarize(HostFitness.mean = mean(HostFitness),
            HostMicrobeFitness.mean = mean(HostMicrobeFitness),
            diversity.mean = mean(diversity, na.rm = TRUE),
            mean_reproductive_output.mean = mean(mean_reproductive_output),
            EnvMicrobeFitness.mean = mean(EnvMicrobeFitness),
            LogFoldChangeMicrobe.mean = mean(LogFoldChangeMicrobe),
            EnvCond.mean = mean(EnvCond),
            LogFoldChangeHost.mean = mean(LogFoldChangeHost),
            AdjRSq.mean = mean(adj_rsq),
            Coeff.mean = mean(host_microbe_coeff, na.rm = TRUE),
            LogFold_ngen1_Host.mean = mean(LogFold_ngen1_Host),
            LogFold_ngen1_Microbe.mean = mean(LogFold_ngen1_Microbe),
            BDiv.mean = mean(host_bdiv, na.rm = TRUE)
  )

all_combinations_evi$VertConOrig<-all_combinations__orig$Var1
all_combinations_evi$NewLookup<-paste(all_combinations_evi$Var1,all_combinations_evi$Var7)
effect_VI_data_means$NewLookup<-paste(effect_VI_data_means$VertCon,effect_VI_data_means$EnvType)
effect_VI_data_means$NewLookup<-paste0(effect_VI_data_means$NewLookup,".1gen")
#effect_VI_data_means$VertCon_200Gen<-all_combinations_evi$VI200Gen[]
effect_VI_data_means$VertCon<-all_combinations_evi$VertConOrig[match(effect_VI_data_means$NewLookup,all_combinations_evi$NewLookup)]
effect_VI_data_means$VertConOrig<-all_combinations_evi$Var1[match(effect_VI_data_means$NewLookup,all_combinations_evi$NewLookup)]
#effect_VI_data_means$VI<-lookupdf$VI[match(effect_VI_data_means$VertCon,round(as.numeric(lookupdf$Val),8))]
#effect_VI_data_means$EffectiveVI<-(as.character(round(as.numeric(multigen_data_means_200gen$VertCon),5)))


#effect_VI_data_means$EffectiveVI<-effect_VI_data_means$VertCon
#effect_VI_data_means$EffectiveVI<-(as.character(round(as.numeric(effect_VI_data_means$EffectiveVI),5)))
#effect_VI_data_means$VertCon<-
effect_VI_data_means$GroupingColumn<-paste(#effect_VI_data_means$EnvType,
                                           effect_VI_data_means$VertCon,
                                           effect_VI_data_means$Gen)
multigen_data_means_200gen<-general_tm_1_10_data_means %>%
  filter(NGens=="200gen") %>% filter(G1=="0")

multigen_data_means_200gen$GroupingColumn<-paste(#multigen_data_means_200gen$EnvType,
                                                 multigen_data_means_200gen$VertCon,
                                                 multigen_data_means_200gen$Gen)
multigen_data_means_200gen$EVI_HostFitness<-effect_VI_data_means$HostFitness.mean[match(multigen_data_means_200gen$GroupingColumn,
                                                                                        effect_VI_data_means$GroupingColumn)]
multigen_data_means_200gen$EVI_MicrobeFitness<-effect_VI_data_means$HostMicrobeFitness.mean[match(multigen_data_means_200gen$GroupingColumn,
                                                                                                  effect_VI_data_means$GroupingColumn)]

#multigen_data_means_200gen$VertCon<-lookupdf$Val[match(multigen_data_means_200gen$VertCon,lookupdf$InitVert)]
#multigen_data_means_200gen$VertCon<-as.character(round(as.numeric(multigen_data_means_200gen$VertCon),8))

multigen_data_means_200gen<-multigen_data_means_200gen[,c("Gen","VertCon","HostMicrobeFitness.mean","HostFitness.mean","EnvType","GroupingColumn")]
multigen_data_means_200gen$Group<-"Multigen"
effect_VI_data_means<-effect_VI_data_means[,c("Gen","VertCon","HostMicrobeFitness.mean","HostFitness.mean","EnvType","GroupingColumn")]
effect_VI_data_means$Group<-"EVI"
effect_VI_data_means$VertCon<-as.character(effect_VI_data_means$VertCon)
effect_VI_data_means<-rbind(effect_VI_data_means,multigen_data_means_200gen)
colnames(effect_VI_data_means)[3]<-"Microbe Fitness"
colnames(effect_VI_data_means)[4]<-"Host Fitness"
effect_VI_data_means$Group<-ifelse(effect_VI_data_means$Group == "EVI","1gen","200gen")
effect_VI_data_means$VertCon<-as.character(round(as.numeric(effect_VI_data_means$VertCon),4))
effect_VI_data_means<-effect_VI_data_means %>% mutate(Group=factor(Group,levels=c("1gen","200gen")))
all_combinations_evi$NewLookup<-paste(all_combinations_evi$Var7,all_combinations_evi$VertConOrig)
effect_VI_data_means$NewLookup<-paste(paste0(effect_VI_data_means$EnvType,".1gen"),effect_VI_data_means$VertCon)
effect_VI_data_means$VertConOrig<-all_combinations_evi$Var1[match(effect_VI_data_means$NewLookup,all_combinations_evi$NewLookup)]
effect_VI_data_means$VertConOrig<-round(effect_VI_data_means$VertConOrig,5)
effect_VI_data_means$VertConOrig<-as.factor(effect_VI_data_means$VertConOrig)
#effect_VI_data_means$VertConEffect<-
efi_hostfit<-effect_VI_data_means %>%
  filter(Gen >= 1500L & Gen <= 1500L) %>%
  ggplot() +
  aes(
    x = Group,
    y = VertCon,
    fill = `Microbe Fitness`
  ) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu", direction = -1) +
  theme_minimal() +
  facet_wrap(vars(EnvType), ncol = 6,
             labeller = labeller(EnvType = c(
               "highac" = "h",
               "highac_incmean" = "i",
               "highac_incvar" = "j",
               "lowac" = "k",
               "lowac_incmean" = "l",
               "lowac_incvar" = "m"
             ))) + xlab(NULL) + 
  ylab("Base Vertical Inheritance") +
  theme(
    strip.text.x = element_text(angle = 0, hjust = 0, face = "bold"), # Bold facet labels
    axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis tick labels
    panel.border = element_blank(),axis.line.y = element_blank(),#(color = "black", fill = NA, size = 1), # Adds full axis lines
    panel.spacing = unit(1, "lines"), # Adds space between panels if necessary
    axis.ticks.length = unit(0.3, "cm"), # Adjust the tick length if necessary
   # axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    strip.background = element_blank()
  )  

efi_microbefit<-effect_VI_data_means %>%
  filter(Gen >= 1500L & Gen <= 1500L) %>%
  ggplot() +
  aes(
    x = Group,
    y = VertCon,
    fill = `Host Fitness`
  ) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu", direction = -1) +
  theme_minimal() +
  facet_wrap(vars(EnvType), ncol = 6,
             labeller = labeller(EnvType = c(
               "highac" = "b",
               "highac_incmean" = "c",
               "highac_incvar" = "d",
               "lowac" = "e",
               "lowac_incmean" = "f",
               "lowac_incvar" = "g"
             )))+ xlab(NULL) + 
  ylab("Base Vertical Inheritance") +
  theme(
    strip.text.x = element_text(angle = 0, hjust = 0, face = "bold"), # Bold facet labels
    axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis tick labels
    panel.border = element_blank(),axis.line.y = element_blank(),#(color = "black", fill = NA, size = 1), # Adds full axis lines
    panel.spacing = unit(1, "lines"), # Adds space between panels if necessary
    axis.ticks.length = unit(0.3, "cm"), # Adjust the tick length if necessary
    #axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    strip.background = element_blank()
  )  

viests$Group<-"x"
efi_relationship<-ggplot(viests) +
  aes(x = Var2, y = value, colour = InitVert) +
  geom_line(linewidth=2) +
  scale_color_hue(direction = 1) +
  theme_classic2() + xlab("Microbial Generations within a host") + 
  ylab("Base Vertical Inheritance")  + 
  facet_wrap(vars(Group), ncol = 1,
             labeller = labeller(Group = c(
               "x" = "a"))) +
  theme(
    strip.text.x = element_text(angle = 0, hjust = 0, face = "bold"), # Bold facet labels
    axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis tick labels
    panel.border = element_blank(),axis.line.y = element_blank(),#(color = "black", fill = NA, size = 1), # Adds full axis lines
    panel.spacing = unit(1, "lines"), # Adds space between panels if necessary
    axis.ticks.length = unit(0.3, "cm"), # Adjust the tick length if necessary
    #axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    strip.background = element_blank()
  )  +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,linewidth=1.8) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,linewidth=1.8) 



library(dplyr)

# Summarize the data to calculate mean and standard deviation
viests_summary <- viests %>%
  group_by(Var2, InitVert) %>%
  summarize(mean_value = mean(value),
            sd_value = sd(value),
            .groups = 'drop') 

efi_relationship <- ggplot(viests_summary) +
  aes(x = Var2, y = mean_value, group = InitVert, colour = InitVert) +
  geom_line(linewidth = 2) + 
  geom_ribbon(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), 
              alpha = 0.2, fill = "grey") + 
  scale_color_hue(direction = 1) +
  theme_classic2() + 
  xlab("Microbial Generations within a host") + 
  ylab("Effective Vertical Inheritance")  + 
  #facet_wrap(vars(InitVert), ncol = 1,
  #           labeller = labeller(InitVert = c(
  #             "x" = "a"))) +
  theme(
    strip.text.x = element_text(angle = 0, hjust = 0, face = "bold"), # Bold facet labels
    axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis tick labels
    panel.border = element_blank(), axis.line.y = element_blank(),
    panel.spacing = unit(1, "lines"), # Adds space between panels if necessary
    axis.ticks.length = unit(0.3, "cm"), # Adjust the tick length if necessary
    strip.background = element_blank()
  )  +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth=1.8) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth=1.8)

efi_relationship


efi_heats<-cowplot::plot_grid(efi_microbefit,efi_hostfit,nrow=2,align = "v")

efi_plot<-cowplot::plot_grid(efi_relationship + theme(legend.position = "none") ,efi_heats,ncol=2,rel_widths = c(2,3))
efi_plot

