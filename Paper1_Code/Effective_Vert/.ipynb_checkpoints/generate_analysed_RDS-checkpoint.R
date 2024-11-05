
library(readr)
library(dplyr)
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
  # print("Issue in extractfitnesses")
  gendata=data$GenData
  data2use=lapply(gendata,function(q){
    env_microbe_fitness_changes<-mean((q$env_fits))#/x[[1]]$GenData[[500]]$env_fits))
    host_microbe_fitness_changes<-mean((q$MicrobeFitness))#/x[[1]]$GenData[[500]]$MicrobeFitness))
    host_fitness_changes<-mean((q$HostFitness))#/x[[1]]$GenData[[500]]$HostFitness))
    gen_of_int_fitness_host<-mean(gendata[[rel_log_gen]]$HostFitness)
    gen_of_int_fitness_microbe<-mean(gendata[[rel_log_gen]]$MicrobeFitness)
    host_microbe_coeff=lm(q$MicrobeFitness~q$HostFitness)
    host_microbe_coeff<-summary(host_microbe_coeff)#$coefficients
    adj_rsq<-host_microbe_coeff$adj.r.squared
    host_microbe_coeff<-as.data.frame(host_microbe_coeff$coefficients)$Estimate[2]
    host_diversity<-mean(q$Diversity)
    host_bdiv<-mean(q$BrayDiv)
    data.frame(EnvMicrobeFitness=env_microbe_fitness_changes,
               HostMicrobeFitness=host_microbe_fitness_changes,
               HostFitness=host_fitness_changes,
               diversity=host_diversity,
               gen_of_int_fitness_microbe=gen_of_int_fitness_microbe,
               gen_of_int_fitness_host=gen_of_int_fitness_host,
               host_microbe_coeff=host_microbe_coeff,
               adj_rsq=adj_rsq,
               host_bdiv=host_bdiv)
  }
  )
  do.call("rbind",data2use)
}

generate_df_genrich<-function(file,env_file,output){
  file<-read_rds(file)
  env_file<-read_rds(env_file)
  env_cond<-unlist(lapply(env_file,function(x){x[,1]}))
  env_df<-data.frame(env_cond,names(env_cond))
  
  #file<-unlist(file,recursive = F)
  gen_rich_varenv<-lapply(file,function(x){
    genetic_richness<-unlist(lapply(x$GenData,function(q){num_unique_values(q$HostMicrobeOptima)}))
    avg_dists<-unlist(lapply(x$GenData,function(q){mean(dist(q$HostMicrobeOptima))}))
    mean_host_prefoptima=unlist(lapply(x$GenData,function(q){mean((q$HostMicrobeOptima))}))
    mean_reproductive_output<-unlist(lapply(x$GenData,function(q){mean((q$microbe_samplingprob))}))
    return(data.frame(AverageDistance=avg_dists,Richness=genetic_richness,mean_host_prefoptima=mean_host_prefoptima,Gen=1:length(mean_host_prefoptima),mean_reproductive_output=mean_reproductive_output))
  })
  gen_rich_varenv<-do.call("rbind",gen_rich_varenv)
  gen_rich_varenv$NGens<-1#stringr::word(rownames(gen_rich),2,sep="[.]")
  fitdata<-lapply(file,function(x){extract_fitnesses(x,200)})
  fitdata<-do.call("rbind",  fitdata)
  gen_rich_varenv<-cbind(gen_rich_varenv,fitdata)
  
  
  #gen_rich_varenv$LookupVal<-stringr::word(rownames(gen_rich_varenv),1,sep=".X")
  
  gen_rich_varenv$LookupVal<-sub(".*?MicrobeSel[^_]*_", "", rownames(gen_rich_varenv))
  gen_rich_varenv$LookupVal <- sub("((\\.[^.]*){1})\\..*", "\\1", gen_rich_varenv$LookupVal)
  gen_rich_varenv$Group_Simulation<-gen_rich_varenv$LookupVal
  
  #gen_rich_varenv$LookupVal <- sub("((\\.[^.]*){1})\\.", "\\1", gen_rich_varenv$LookupVal)
  #gen_rich_varenv$LookupVal<-sub('.', '', gen_rich_varenv$LookupVal)
  
  gen_rich_varenv$LookupVal<-paste0(gen_rich_varenv$LookupVal,gen_rich_varenv$Gen)
  
  gen_rich_varenv$EnvCond<-env_df$env_cond[match(gen_rich_varenv$LookupVal,env_df$names.env_cond.)]
  #gen_rich_varenv$Group_Simulation<-unlist(group_sim)#c(rep("env_list_increasing_var",12000),rep("increasing_env",12000),rep("Static",12000))
  gen_rich_varenv$EnvType<-stringr::word(gen_rich_varenv$Group_Simulation,1,sep="[.]")
  gen_rich_varenv$NGens<-stringr::word(gen_rich_varenv$Group_Simulation,2,sep="[.]")
  gen_rich_varenv$VertCon<-stringr::word(rownames(gen_rich_varenv),2,sep="X")
  gen_rich_varenv$VertCon<-stringr::word(gen_rich_varenv$VertCon,1,sep="_Y")
  
  
  gen_rich_varenv$Host2Env<-stringr::word(rownames(gen_rich_varenv),5,sep="_")
  gen_rich_varenv$Host2Env<-stringr::str_sub(gen_rich_varenv$Host2Env,-1,-1)
  gen_rich_varenv$Env2Env<-stringr::word(rownames(gen_rich_varenv),6,sep="_")
  gen_rich_varenv$Env2Env<-stringr::str_sub(gen_rich_varenv$Env2Env,-1,-1)
  
   gen_rich_varenv<-gen_rich_varenv %>%
    mutate(NGens=factor(NGens,levels=c("1gen","2gen","5gen","10gen","20gen","50gen","100gen","200gen")))
  
  gen_rich_varenv$LogFoldChangeMicrobe=log(gen_rich_varenv$HostMicrobeFitness/gen_rich_varenv$gen_of_int_fitness_microbe)
  gen_rich_varenv$LogFoldChangeHost=log(gen_rich_varenv$HostFitness/gen_rich_varenv$gen_of_int_fitness_host)
  gen_rich_varenv$MicrobeSel<-sub(".*?([^_.]*MicrobeSel[^_.]*).*", "\\1", rownames(gen_rich_varenv)[grep("MicrobeSel", rownames(gen_rich_varenv))])
  gen_rich_varenv$HostSel<-sub(".*?([^_.]*HostSel[^_.]*).*", "\\1", rownames(gen_rich_varenv)[grep("HostSel", rownames(gen_rich_varenv))])
  readr::write_rds(gen_rich_varenv,output)
}



args <- commandArgs(trailingOnly = TRUE)

# Check if three arguments are provided
if (length(args) != 3) {
  stop("Please provide three arguments: file, envdata, and outputfile.")
}

# Assign arguments to variables
data_file_list <- args[1]
env_file_list <- args[2]
new_filename <- args[3]
process_element <- function(i) {
  generate_df_genrich(data_file_list, env_file_list, new_filename)
}

# Use mclapply to apply the function in parallel
library(parallel)
replicate_data <- mclapply(1:length(data_file_list), process_element, mc.cores = 1)